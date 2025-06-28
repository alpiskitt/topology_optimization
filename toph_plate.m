%%%% Laplace heat topology optimization %%%
function toph_plate(nelx,nely,volfrac,penal,rmin,ft,method)
% INITIALIZE
x(1:nely,1:nelx) = volfrac; 
loop = 0; 
change = 1.0;
%%%%% --- BUILD FILTER --------------------------------------------------- %%
%  iH, jH   index vectors, sH = weight,  Hs = row sums  (just like Sigmund’01)
iH = [];  jH = [];  sH = [];
for i = 1:nelx
  for j = 1:nely
    e1 = (i-1)*nely + j;
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        e2 = (k-1)*nely + l;
        fac = rmin - sqrt((i-k)^2 + (j-l)^2);
        if fac > 0
          iH(end+1) = e1;
          jH(end+1) = e2;
          sH(end+1) = fac;
        end
      end
    end
  end
end
H  = sparse(iH,jH,sH);
Hs = sum(H,2);                      % column vector, length = nelx*nely
% --- MMA constants ----------------------------------------------------- %
m   = 1;                           % number of constraints (volume only)
n   = nelx*nely;                   % number of design variables
xmin  = zeros(n,1);  xmax = ones(n,1);
low   = [];  upp = [];             % moving asymptotes (first call = auto)
xval  = volfrac*ones(n,1);         % MMA column design vector
xold1 = xval;   xold2 = xval;      % two previous designs
a0 = 1;  a = zeros(m,1);
c_mma = 5e6*ones(m,1);   % << instead of 1e3
d     = zeros(m,1);      % leave quadratic term zero


% ---------------------------------------------------------------
% nested helper: run FE-analysis and return compliance for a given density field
function comp = compliance(xPhys)
    sK  = reshape(KE(:)*(1e-6 + xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K   = sparse(iK,jK,sK);  K = (K+K')/2;
    U   = zeros(2*(nely+1)*(nelx+1),1);
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    ce  = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    comp = sum(sum((Emin + xPhys.^penal*(E0-Emin)).*ce));
end
% pick one design-variable index to test (here dead-centre element)
testEl = sub2ind([nely nelx], round(nely/2), round(nelx/2));
%% ----  MMA bookkeeping -------------------------------------------------
m     = 1;                     % number of constraints  (here: volume <= V*) only
n     = nelx*nely;             % number of design variables
xmin  = zeros(n,1);            % lower bounds on x
xmax  = 1*ones(n,1);           % upper bounds
low   = [];                  % initialise asymptotes
upp   = [];
xval  = volfrac*ones(n,1);     % column vector version of x
xold1 = xval;                  % two previous designs (zero-change start)
xold2 = xval;
a0 = 1;  a = zeros(m,1);       % standard MMA constants
c  = 1000*ones(m,1);
d  = zeros(m,1);
%%
% START ITERATION
while change > 0.01  
  loop = loop + 1;
  xold = x;
  if ft == 1
      xPhys = x;
  else
      xPhys = reshape((H*x(:))./Hs , nely, nelx);
  end

% FE-ANALYSIS
  L   = 0.10;                       % plate size  [m]
  a   = 0.015;                      % half die-width
  q0  = 1e6;                        % heat flux   [W/m²]
  k_s = 1;                          % normalise solid conductivity to 1
  [~,~,U] = FE(nelx,nely,xPhys,penal,q0,a,L,k_s);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;
  c = 0.;
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      Ue = U([n1; n2; n2+1; n1+1],1);
      val = xPhys(ely,elx);                    
      kappa = 0.001 + 0.999*val^penal;      % same factor
      c  = c + kappa * Ue' * KE * Ue;
      dc(ely,elx) = -0.999*penal*val^(penal-1) * Ue' * KE * Ue;
    end
  end
%%%%% --- FILTERING OF SENSITIVITIES ----------------------------------- %%
if ft == 1                      % traditional “sensitivity” filter
    dc(:) = (H*(x(:).*dc(:)))./Hs ./ max(1e-3,x(:));
else                            % density filter (ft == 2)
    dc(:) = H*dc(:)./Hs;        % already satisfies chain rule
end
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
dv = ones(nely,nelx);
% ---- FILTERED sensitivities ready:  dc ,  dv -------------------------- %
% ---- volume derivative -----------------------------------------------
if ft == 2                      % density filter: dv must be filtered too
    dv(:) = H*(dv(:)./Hs); 
end
switch method
  %======================================================================
  case "oc"
      % --- standard optimality-criteria update ------------------------ %
      l1 = 0;  l2 = 1e5;  move = 0.2;
      while (l2-l1) > 1e-4
          lmid = 0.5*(l2+l1);
          xnew = max(0.001, ...
                max(x-move, ...
                min(1.0, min(x+move, x.*sqrt(-dc./dv/lmid)))));
          if ft==1
              xPhys = xnew;
          else
              xPhys = reshape((H*xnew(:))./Hs , nely, nelx);
          end
          if mean(xPhys(:)) - volfrac > 0,  l1 = lmid;  else  l2 = lmid; end
      end

  %======================================================================
  case "mma"
      % --- objective & gradient in column-vector form ----------------- %
  
      f0val = c;            df0dx = dc(:);
      % single constraint: g(x) = mean(xPhys) - volfrac <= 0
      %fval  = (mean(xPhys(:)) - volfrac)/volfrac;

      % ----- volume constraint ----------------------------------------------
      Vtarget = volfrac*nelx*nely;            % 480 for a 40×40 mesh
      g       = sum(xPhys(:)) - Vtarget;      % f_i(x)  ( >0  ⇒ violation)

    if ft==1
        dgdX =  ones(n,1);                  % ∂g/∂x   (column)
    else
        dgdX =  dv(:);                      % density-filter derivative
    end

    % --- make its gradient comparable to the objective --------------------
    scale =  max(abs(df0dx)) / max(abs(dgdX));   % O(100–1000)
    g      =  g      * scale;
    dgdX   =  dgdX   * scale;
    
    fval   =  g;           % MMA wants row-vector(s) later
    dfdx   =  dgdX.';      % (1×n)

      [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
          mmasub(m,n,loop,x(:),xmin,xmax,xold1,xold2, ...
                  f0val,df0dx,fval,dfdx,low,upp,a0,a,c_mma,d);
      move = 0.02;
      xmma = max(x(:)-move, min(x(:)+move, xmma));  % clip

    % ---- store *clipped* design for next sub-problem -------------------
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;                                 % <-- AFTER clipping
    xnew  = reshape(xval,nely,nelx);
      if ft==1
          xPhys = xnew;
      else
          xPhys = reshape((H*xnew(:))./Hs , nely, nelx);
      end
  %======================================================================
     case "fdcheck"
          % we need to do finite difference
          if loop == 1
            c0    = c;              % <-- call it c0 here …
            g_adj = dc(testEl);
            xkeep = x;
            relerr_log = [];  dp_log = [];
            kFD = 0;
          end
          kFD = kFD +1 ;
    
    % ---- choose perturbation size (10^{-loop}) -------------------
          delta_p = 10^(-kFD); 
    
    % ---- perturb ONE variable, keep inside [0,1] -----------------
         xpert        = xkeep;           % start from original design each time
         xpert(testEl)= min(1-1e-9, xpert(testEl)+delta_p);
        
    % ---- filter to physical densities if needed ------------------
         if ft==1
            xPhys_pert = xpert;
         
         else
            tmp          = (H*xpert(:))./Hs;        % density filter result (vector)
            xPhys_pert   = reshape(tmp, nely, nelx); % back to nely-by-nelx matrix
         end         
                
    % ---- compliance with perturbed design ------------------------
         c_pert = compliance(xPhys_pert);
    
    % ---- finite-difference estimate and relative error -----------
         g_fd   = (c_pert - c0)/delta_p;
         relerr = abs(g_fd - g_adj) / max(1e-30, abs(g_adj));
    
    % ---- store and print -----------------------------------------
         relerr_log(end+1) = relerr;
         dp_log   (end+1)  = delta_p;
         fprintf('FD-check  Δ=%-10.1e   grad_fd=%-10.3e   err=%-10.3e\n',...
                delta_p, g_fd, relerr);
    
    % ---- stop when Δ below machine precision or 6 decades tested -
        if kFD >= 6             % 10^-1 … 10^-6 done
            figure, loglog(dp_log, relerr_log, '-o');
            xlabel('\delta_p'); ylabel('relative error');
            title('FD sensitivity check');
            grid on, axis tight
            break                % leave the WHILE completely
        end
         
         % no design update (x stays xkeep), just continue to next Δ
         change = 1;
         xnew   = xkeep;

        continue
  otherwise
      error("method must be ""oc"" or ""mma""")
end

% -------- convergence & bookkeeping ----------------------------------- %
change = max(abs(xnew(:)-x(:)));
x      = xnew;
% ---- diagnostics ------------------------------------------------------
fprintf(' It.:%4i  Obj.:%10.4f  Vol.:%6.3f  ch.:%6.3f\n', ...
        loop, c, mean(xPhys(:)), change);

colormap(gray);
imagesc(1 - xPhys);         % invert colours: black = material
axis equal off;  drawnow;

end 
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FE-ANALYSIS  (heat conduction) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,F,U] = FE(nelx,nely,xPhys,penal,q0,a,L,k_s)
    [KE] = lk;                            % 4-node Laplace element
    K = sparse((nelx+1)*(nely+1), (nelx+1)*(nely+1));
    F = sparse((nelx+1)*(nely+1),1);      % right-hand side = heat sources
    h = L/nelx;                           % element size (square mesh)

    % loop over elements --------------------------------------------------
    for elx = 1:nelx
        for ely = 1:nely
            n1 = (nely+1)*(elx-1)+ely;
            n2 = (nely+1)* elx   +ely;
            edof = [n1; n2; n2+1; n1+1];

            kappa = (0.001 + 0.999*xPhys(ely,elx)^penal);   % interpolated k/k_s
            K(edof,edof) = K(edof,edof) + kappa * KE;

            % --- heat generation in the central die ---------------------
            xc = (elx-0.5)*h - L/2;      % element centre (x,y)
            yc = (ely-0.5)*h - L/2;
            if abs(xc) <= a && abs(yc) <= a
                F(edof) = F(edof) + q0 * h^2 / 4;   % 4 nodes share q·area
            end
        end
    end

    % --- Dirichlet boundary (outer frame at T=0) -------------------------
    frameNodes = unique([1:(nely+1), ...                     % left
                         (nely+1): (nely+1): (nelx)*(nely+1)+1, ... % top
                         (nelx)*(nely+1)+1 : (nelx+1)*(nely+1), ...  % right
                         1:(nely+1):(nelx)*(nely+1)+1 ]);    % bottom
    alldofs = 1:(nelx+1)*(nely+1);
    freedofs = setdiff(alldofs,frameNodes);

    % --- solve -----------------------------------------------------------
    U = zeros((nelx+1)*(nely+1),1);
    U(freedofs) = K(freedofs,freedofs) \ F(freedofs);
    % fixed nodes already at zero
end

%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
KE = [ 2/3 -1/6 -1/3 -1/6
      -1/6 2/3 -1/6 -1/3
      -1/3 -1/6 2/3 -1/6
      -1/6 -1/3 -1/6 2/3];
end
end