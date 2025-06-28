%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function [x, c] = top88_pbm8(nelx,nely,volfrac,penal,rmin,ft,method,omega,x0)
close all
%% CHOICE
passive = "no";
active = "no";
delta_p = 10^-1;
rhomin = 10^-8;

%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-4;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

% The element mass matrix
m0 = [4/9 0 2/9 0 1/9 0 2/9 0 
   0 4/9 0 2/9 0 1/9 0 2/9 
   2/9 0 4/9 0 2/9 0 1/9 0 
   0 2/9 0 4/9 0 2/9 0 1/9 
   1/9 0 2/9 0 4/9 0 2/9 0 
   0 1/9 0 2/9 0 4/9 0 2/9 
   2/9 0 1/9 0 2/9 0 4/9 0 
   0 2/9 0 1/9 0 2/9 0 4/9]/4/(nelx*nely); 

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
% PlotLoadsAndSupports2021(nelx, nely, fixeddofs, F(:,1), 0, 'no', 'no')
% pause

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
if size(x0) == [nely, nelx]
    x = x0;
    if ft == 1
        xPhys = x;
    elseif ft == 2
        xPhys(:) = (H*x(:))./Hs;
    end
else
    x = repmat(volfrac,nely,nelx);
end
xmin = repmat(0,nely,nelx);
xmax = repmat(1,nely,nelx);

% Passive zone
if passive == "yes"
    x_passive = ceil((nelx)/3) + 1 : nelx - ceil((nelx)/3);
    y_passive = ceil((nely)/3) + 1 : nely - ceil((nely)/3);
    
    counter = 1;
    for i = 1:length(x_passive)
        for j = 1:length(y_passive)
            xmax(y_passive(j), x_passive(i)) = 0.05;
            counter = counter + 1;
        end
    end
end

% Active zone
if active == "yes"
    P1 = [1, 1];
    P2 = [nely, 1];
    P3 = [1, nelx];
    
    value_active = 0.9;
    xmin(P1) = value_active;
    xmin(P2) = value_active;
    xmin(P3) = value_active;
end

xPhys = x;
loop = 0;
change = 1;
low = zeros(prod(size(x)), 1);
upp = zeros(prod(size(x)), 1);

xold2 = x;
xold1 = x;

%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin + xPhys(:)'.^penal*(E0 - Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  mK = reshape(m0(:)*(rhomin + (1 - rhomin)*xPhys(:)'),64*nelx*nely,1);
  M = sparse(iK,jK,mK); M = (M+M')/2;
  S = K - omega^2 * M;
  U(freedofs,1) = S(freedofs,freedofs) \ F(freedofs);
  %% OBJECTIVE FUNCTION, CONSTRAINTS, SENSITIVITY ANALYSIS
  % Objective
  UU = U(:,1);
  ceKU = reshape(sum((UU(edofMat)*KE).*UU(edofMat),2),nely,nelx);
  ceMU = reshape(sum((UU(edofMat)*m0).*UU(edofMat),2),nely,nelx);
  c = norm(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ceKU)) - ...
      omega^2*sum(sum(((rhomin + (1 - rhomin)*xPhys).*ceMU))));
  
  f0val = c;
  % Gradient of objecti
  % Adjoint
  gamma = 0.5 * sign(conj(U(:,1))'*F(:));
  U(freedofs,2) = -gamma*U(freedofs,1);
  UL = U(:,2);
  ceK = reshape(sum((UL(edofMat)*KE).*UU(edofMat),2),nely,nelx);
  ceM = reshape(sum((UL(edofMat)*m0).*UU(edofMat),2),nely,nelx);
  
  dce = (penal*(E0-Emin)*xPhys.^(penal-1)).*ceK - omega^2*(1-rhomin).*ceM;
  dc = 2* real(dce);
  
  % Constraints
  vol = mean(xPhys(:))/volfrac - 1.0;
  fval = vol;

  % Gradient of constraints
  dv = (1/(volfrac * nelx * nely)) * ones(nely,nelx);

  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end

  df0dx = dc;

  dfdx = dv;
  switch method
      case "oc"  
            %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
            l1 = 0; l2 = 1e9; move = 0.2;
            while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1);
            xnew = max(xmin,max(x-move,min(xmax,min(x+move,x.*sqrt(-dc./dv/lmid)))));
            if ft == 1
              xPhys = xnew;
            elseif ft == 2
              xPhys(:) = (H*xnew(:))./Hs;
            end
            if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
            end
            change = max(abs(xnew(:)-x(:)));
            x = xnew;
      case "mma"
          %% MMA UPDATE
          m_mma = 1;
          n_mma = nelx*nely;
          a0_mma = 1;
          a_mma = zeros(m_mma,1);
          c_mma = 1000 * ones(m_mma,1);
          d_mma = 0;

          [xnew, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = mmasub(m_mma, n_mma, loop, ...
                    vec2mat(x), vec2mat(xmin), vec2mat(xmax),...
                    vec2mat(xold1), vec2mat(xold2), ...
                f0val, vec2mat(df0dx), fval, vec2mat(dfdx)', low, upp, ...
                a0_mma, a_mma, c_mma, d_mma);
          xnew = mat2vec(xnew, nely, nelx);

          change = max(abs(xnew(:)-x(:)));

          xold2 = xold1;
          xold1 = x;
          x = xnew;
          if ft == 1
            xPhys = xnew;
          elseif ft == 2
            xPhys(:) = (H*x(:))./Hs;
          end
      case "fdcheck"
          for i = 1:nely
              for j = 1:nelx
                xmod = x;
                xmod(i, j) = xmod(i, j) + delta_p;
                xmodPhys = xmod;
                if ft == 1
                 xmodPhys = xmod;
                elseif ft == 2
                 xmodPhys(:) = (H*xmod(:))./Hs;
                end
                %% FE-ANALYSIS
                  sK = reshape(KE(:)*(Emin + xmodPhys(:)'.^penal*(E0 - Emin)),64*nelx*nely,1);
                  K = sparse(iK,jK,sK); K = (K+K')/2;
                  mK = reshape(m0(:)*(rhomin + (1 - rhomin)*xmodPhys(:)'),64*nelx*nely,1);
                  M = sparse(iK,jK,mK); M = (M+M')/2;
                  S = K - omega^2 * M;
                  U(freedofs,1) = S(freedofs,freedofs) \ F(freedofs);
                  %% OBJECTIVE FUNCTION, CONSTRAINTS, SENSITIVITY ANALYSIS
                  % Objective
                  UU = U(:,1);
                  ceKU = reshape(sum((UU(edofMat)*KE).*UU(edofMat),2),nely,nelx);
                  ceMU = reshape(sum((UU(edofMat)*m0).*UU(edofMat),2),nely,nelx);
                  cmod = norm(sum(sum((Emin+xmodPhys.^penal*(E0-Emin)).*ceKU)) - ...
                      sum(sum((omega^2*(rhomin + (1 - rhomin)*xmodPhys).*ceMU))));
  
                  % Gradient of objective
                  dcmod(i,j) = (cmod - c) / delta_p;
              end
          end
          err(loop) = norm(vec2mat(dcmod) - vec2mat(dc)) / norm(vec2mat(dc));
          delta_p = delta_p / 10;
          if loop > 10
              change = 0.005;
          end
      otherwise
          error("method not implemented")
          
  end

  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; 
  title(strcat("\omega = ", num2str(omega)), "interpreter", "tex");
  drawnow;

end
switch method
    case "fdcheck"
        loglog(10.^[-1:-1:-length(err)], err, '-o');
    otherwise
end

end

function v = vec2mat(A)
    v = reshape(A, [prod(size(A)), 1]);
end

function A = mat2vec(v, nrow, ncol)
    if nrow*ncol ~= length(v)
        error("size not compatible");
    else
        A = reshape(v, [nrow, ncol]);
    end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%