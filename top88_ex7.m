%% MATERIAL PROPERTIES
function top88_ex7(method)
nelx = 300;
nely = 100;
volfrac = 0.3;
penal = 3.0;
rmin = 1.3;
ft = 2;
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
%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
first_node = 1;
top_right = 2*nelx*(nely+1)+1;

F = sparse([first_node top_right],[1 2],[1 -1],2*(nely+1)*(nelx+1),2); 
U = zeros(2*(nely+1)*(nelx+1),2);
% --- boundary conditions -------------------------------------------------
% 1) uy = 0 along the whole left edge  (roller support)
fixeddofs = 2 : 2*(nely+1) : 2*(nely+1)*(nelx+1);

% 2) clamp ONLY the bottom-left node  ──► change happens here
cornerNode          = nely + 1;                     % node number at bottom-left
fixeddofs_corner    = [ 2*cornerNode-1 , 2*cornerNode ];   % = [ux  uy] of that node
% If you want a roller instead of a clamp, keep only the x-dof:
% fixeddofs_corner  = 2*cornerNode-1;

% merge the two sets
fixeddofs           = union(fixeddofs , fixeddofs_corner);

% assemble the free-dof list
alldofs  = 1 : 2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs , fixeddofs);
%% VIEW NODES AND SUPPORTS
% F = sparse([first_node top_right],[1 1],[1 -1],2*(nely+1)*(nelx+1),2); 
% PlotLoadsAndSupports2021(nelx, nely, fixeddofs, F, [], 'yes', 'yes')
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

%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x;
loop = 0;
change = 1;
rel_errors = zeros(1,10); %fdcheck
delta_ps = zeros(1,10);   %fdcheck
ele_pert = nely*4;     %fdcheck/must be between 1 and nelx*nely
%% PASSIVE ELEMENTS
passive = zeros(nely,nelx);  
passive(nely,1) = 2; % material at constraint
passive(1,1) = 2;
passive(1,nelx) = 2;
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

%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  K(first_node, first_node) = K(first_node, first_node) + 0.01;
  K(top_right, top_right) = K(top_right, top_right) + 0.01;
  U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:); % we need to solve for both lambda and u
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  U1 = U(:,1);  
  U2 = U(:,2);  
  ce = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2), nely,nelx);  
  C = U(top_right, 1); 
  dc = penal*(E0-Emin)  *xPhys.^(penal-1).*ce;  
  dv = ones(nely,nelx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
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
            arg  = max(1e-9 , -dc./dv/lmid);   % keep it positive
            xnew = max(0.001, max(x-move, min(1, min(x+move, x .* arg.^0.3))));
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
          % ----------------------------------------------------------
          % Finite-difference gradient test … plus a log-log plot
          % ----------------------------------------------------------
          ndiff        = 10;                      % # elements to sample
          elems2chk    = randperm(nelx*nely,ndiff);
          dp_vec       = logspace(-1,-7,7);       % perturbation sizes
          relerr_vec   = zeros(size(dp_vec));     % store mean |ε| for each Δp

          for idp = 1:numel(dp_vec)
              delta_p = dp_vec(idp);              % this step size
              relerr  = zeros(1,ndiff);           % element-wise errors

              for k = 1:ndiff
                  e  = elems2chk(k);

                  % -- build perturbed design ------------------------
                  x_fd        = x(:);
                  x_fd(e)     = min(1 , x_fd(e)+delta_p);

                  if ft == 1
                      xPhys_fd = reshape(x_fd , nely , nelx);
                  else
                      xPhys_fd = reshape((H*x_fd)./Hs , nely , nelx);
                  end

                  % -- FEA for perturbed design ----------------------
                  sK_fd = reshape(KE(:).*(Emin + xPhys_fd(:)'.^penal*(E0-Emin)), ...
                                   64*nelx*nely,1);
                  K_fd  = sparse(iK,jK,sK_fd);  K_fd = (K_fd+K_fd')/2;
                  K_fd(first_node,first_node) = K_fd(first_node,first_node) + 0.01;
                  K_fd(top_right ,top_right)  = K_fd(top_right ,top_right)  + 0.01;

                  U_fd            = zeros(size(U));
                  U_fd(freedofs,:)= K_fd(freedofs,freedofs)\F(freedofs,:);
                  C_fd            = U_fd(top_right,1);

                  % -- FD derivative & relative error ---------------
                  fd_deriv  = (C_fd - C)/delta_p;
                  ana_deriv = dc(e);
                  relerr(k) = abs(fd_deriv - ana_deriv) / max(1e-12,abs(ana_deriv));
              end

              relerr_vec(idp) = mean(relerr);   % one number per Δp
              fprintf('Δp=%8.1e   mean rel.err = %.3e\n', ...
                       delta_p, relerr_vec(idp));
          end

          % ---------- plot error vs. perturbation size ---------------
          figure, loglog(dp_vec, relerr_vec, '-o');
          xlabel('\delta_p'); ylabel('relative error');
          title('FD sensitivity check');
          grid on, axis tight

          break                      % leave the WHILE completely
          
          
        otherwise
            error("Method not implemented");
    end


  
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
     %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,C, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
end
