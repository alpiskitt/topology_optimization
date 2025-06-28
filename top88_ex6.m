%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function top88_ex6(nelx,nely,volfrac,penal,rmin,ft, method)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
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

% ---------------------------------------------------------------
% nested helper: run FE-analysis and return compliance for a given density field
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
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x;
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
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
          f0val  = c;                         % compliance from above
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
            error("Method not implemented");
    end


  
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end