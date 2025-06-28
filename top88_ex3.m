%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function top88_ex3(nelx,nely,volfrac,penal,rmin,ft)
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
%% DEFINE LOADS AND SUPPORTS (“FULL” BRIDGE WITH POINT OR DISTRIBUTED LOAD)
ndof   = 2*(nely+1)*(nelx+1);   % total number of DOFs (helps readability)

F1 = sparse(ndof,1);            
F2 = sparse(ndof,1);            

U1 = zeros(ndof,1);             
U2 = zeros(ndof,1);             

% --- bottom supports (last row of nodenrs) ---
nL = nodenrs(nely+1,1);
nR = nodenrs(nely+1,nelx+1);
fixeddofs = [2*nL-1 2*nL 2*nR-1 2*nR];

alldofs  = 1:ndof;                %  ← restore this line
freedofs = setdiff(alldofs,fixeddofs);   %  ← …and this one
   
% point loads at 1/3 & 2/3 along the *top* edge (= first row)
i1 = round(nelx/3);    
i2 = round(2*nelx/3);
n1 = nodenrs(1,i1+1); 
n2 = nodenrs(1,i2+1);

F1(2*n1) = -1;                  % ↓ use the new, full-length vectors
F2(2*n2) = -1;
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
%PlotLoadsAndSupports2021(nelx, nely, fixeddofs, F1, [ ], 'yes' , 'yes');
% pause;
%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK);  K = (K+K')/2;

  U1(freedofs) = K(freedofs,freedofs)\F1(freedofs);
  U2(freedofs) = K(freedofs,freedofs)\F2(freedofs);
  ce1 = reshape(sum((U1(edofMat)*KE).*U1(edofMat),2),nely,nelx);
  ce2 = reshape(sum((U2(edofMat)*KE).*U2(edofMat),2),nely,nelx);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  c1  = sum(sum((Emin + xPhys.^penal * (E0-Emin)).*ce1));
  c2  = sum(sum((Emin + xPhys.^penal * (E0-Emin)).*ce2));
  c   = c1 + c2;              % overall objective  (← new)

  dc1 = -penal*(E0-Emin)*xPhys.^(penal-1).*ce1;
  dc2 = -penal*(E0-Emin)*xPhys.^(penal-1).*ce2;
  sumdc = dc1 + dc2;          % combined sensitivity (already in your code)

  dv = ones(nely,nelx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1 %sensitivity filter
    dc1(:) = H*(x(:).*dc1(:))./Hs./max(1e-3,x(:));
    dc2(:) = H*(x(:).*dc2(:))./Hs./max(1e-3,x(:));

  elseif ft == 2  %density filter
    dc1(:) = H*(dc1(:)./Hs);
    dc2(:) = H*(dc2(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end

  sumdc = dc1 + dc2;
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-sumdc./dv/lmid)))));
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    end
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
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

