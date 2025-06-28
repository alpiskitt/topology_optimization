%% The time harmonic dynamic compliance case with MMA and 
function top88_ex8(method)
E0 = 1;
nu = 0.3;
nelx= 300;
nely = 100;
volfrac = 0.5;
penal = 3.0;
rmin = 1.4;
ft = 2;
omega = 0.00;
Emin = 1e-4;
rho0    = 1;                 % solid density                   %%% <<< P8
rhomin  = 1e-8; 
nDof =  2*(nelx+1)*(nely+1);

S      = load('xPhys_static.mat');  % S.xPhys now holds the data
signal = S.xPhys;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nel = nelx*nely;
m0 = (1/36)* [4 0 2 0 1 0 2 0
              0 4 0 2 0 1 0 2
              2 0 4 0 2 0 1 0
              0 2 0 4 0 2 0 1
              1 0 2 0 4 0 2 0
              0 1 0 2 0 4 0 2
              2 0 1 0 2 0 4 0
              0 2 0 1 0 2 0 4];


nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),2);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

function comp = compliance(xPhys)
    sK  = reshape(KE(:)*(Emin + xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
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
n     = nelx*nely;             % number of design variables basically dof
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
%% INITIALIZE ITERATION
%PlotLoadsAndSupports2021(nelx, nely, fixeddofs, F, [], 'yes', 'yes')

if omega == 0
    x = repmat(volfrac,nely,nelx);
    xPhys = x;
else 
    x = signal;
    xPhys = signal;
end

loop = 0;
change = 1;
%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  %% — now build the *global* mass matrix
  % 1) local 8×8 element mass matrix from Appendix B (already in your m0)
  %    note you defined m0 above as an 8×8 and divided by (4*nel)
  rho_e = rhomin + (rho0 - rhomin) * xPhys(:);      % linear mass interpolation should be implemented
  % 3) expand to all 64 entries per element:
  %    make a (nelx*nely × 64) array where each row is m0(:)' times rho_e(e)
  Me_all = rho_e * (m0(:)');                        % size = [nelx*nely, 64]
  sM     = reshape( Me_all.', 64*nelx*nely, 1 );    % vectorized for sparse
  % 4) assemble into the big nDof×nDof sparse matrix
  M = sparse( iK, jK, sM, nDof, nDof );  M = (M + M') * 0.5;     
  
  %% — form the *dynamic stiffness* and solve
  S = K - omega^2 * M;      
  U(freedofs,1) = S(freedofs,freedofs) \ F(freedofs);
  UU = U(:,1);
  ce_K = reshape(sum((UU(edofMat)*KE).*UU(edofMat),2),nely,nelx);
  ce_M = reshape(sum((UU(edofMat)*m0).*UU(edofMat),2),nely,nelx);
  %C = norm(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce_M)) -omega^2* sum(sum(((rhomin+(1 - rhomin)*xPhys).*ce_M))));
  C     = norm(F'*conj(U));
  % Adjoint
  gamma = 0.5*sign(conj(U(:,1))'*F(:));
  U(freedofs,2) = -gamma*U(freedofs,1);
  UA = U(:,2);

  ce_k = reshape(sum((UA(edofMat)*KE).*UU(edofMat),2),nely,nelx);
  ce_m = reshape(sum((UA(edofMat)*m0).*UU(edofMat),2),nely,nelx);

  dce = (penal*(E0-Emin)*xPhys.^(penal-1)).*ce_k - omega^2*(1-rhomin).*ce_m;

  dc = 2*real(dce);

  dv = ones(nely,nelx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1 %sensitivity filter
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2  %density filter
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  switch method 
      case "oc"
        l1 = 0; l2 = 1e9; move = 0.2;
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1);
            xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
            if ft == 1
                xPhys = xnew;
            elseif ft == 2
                xPhys(:) = (H*xnew(:))./Hs;
            end
            if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
        end

      case "mma"
          % objective and gradient in column-vector form
          f0val  = C;                         % compliance from above
          df0dx  = dc(:);                     % sensitivities (dC/dx)
        
          % single volume constraint:  g(x) = mean(xPhys) - volfrac <= 0
          fval   = (mean(xPhys(:)) - volfrac)/volfrac;
          dfdx   = (1/(nelx*nely))*ones(1,n)/volfrac; % dg/dx  (row vector is fine)
        
          [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
                mmasub(m,n,loop,x(:),xmin,xmax,xold1,xold2,   ...
                        f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
        
          xold2 = xold1;
          xold1 = xval;
          xval  = xmma;               % new design from MMA
          xnew  = reshape(xval,nely,nelx);
        
          % physical densities with density filter (unchanged)
          if ft == 1
              xPhys = xnew;
          else
              xPhys(:) = (H*xnew(:))./Hs;
          end
      case "fdcheck"
          % we need to do finite difference
          if loop == 1
            c0    = C;              % <-- call it c0 here …
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
            error("Method not implemented");
    end


  
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,C, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  %colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
  colormap(gray); axis equal; axis tight; axis off;
  Faces = ((edofMat(:,[1 3 5 7])-1)/2)+1;
  [XGrid,YGrid]=meshgrid(0:nelx,nely:-1:0);
  Grid = [XGrid(:),YGrid(:)];
  Deform = [U(1:2:end,1),U(2:2:end,1)];
  patch('Faces', Faces,'Vertices', Grid+Deform*0.05,'FaceColor','Flat','FaceVertexCData',1-x(:),'EdgeColor','none');
  drawnow; clf;

end
%%  FREQUENCY-RESPONSE PLOT  (|uᵀf| vs ω  on log-y scale)
%% ==============================================================

% --- 1.  choose a frequency range ------------------------------
%     an easy way is to span 0 … 2·ω₁ where ω₁ is the 1st
%     (undamped) natural circular frequency of the current
%     structure:

[V,D] = eigs(K(freedofs,freedofs), M(freedofs,freedofs), 1, 'SM');
w1    = sqrt(D(1,1));                     % first eigen-frequency  (rad/s)

nPts  = 400;                              % resolution of the sweep
wVec  = linspace(0, 2*w1, nPts);          % 0 … 2·ω₁
PhiVec = zeros(size(wVec));               % pre-allocate results

% --- 2.  frequency sweep ---------------------------------------
for k = 1:nPts
    w  = wVec(k);
    S  = K - w^2 * M;                     % dynamic stiffness
    U(freedofs) = S(freedofs,freedofs) \ F(freedofs);
    PhiVec(k)   = abs(F.'*U);             % |uᵀf|
end

% --- 3.  plot ---------------------------------------------------
figure('Name','Frequency response');
semilogy(wVec, PhiVec,'LineWidth',1.2);   % log scale in y
xlabel('\omega  [rad/s]');
ylabel('|u^T f|  (dynamic compliance)');
title('Frequency response of current design');
grid on;

% save xPhys_static xPhys
end