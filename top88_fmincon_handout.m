%%%% Compliance topology optimization using FMINCON/L-BFGS
% and general solver comparison
% 
% Please note that the IPOPT optimizer is available as a module 
% https://se.mathworks.com/matlabcentral/fileexchange/53040-ebertolazzi-mexipopt
%
% Density filter is employed for all computations.
%
%%%% Based on:
%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
%
% Version modified for 41591 Spring 2021
function top88_fmincon_handout(nelx,nely,volfrac,penal,rmin)
close all
% Choose the optimizer to run; ipopt, fmincon, mma, oc
optimizer='fmincon';
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
%% STORE FILTER AS A SINGLE SPARSE MATRIX
HsH = spdiags(1.0./Hs,0,nelx*nely,nelx*nely)*H;
%% PACK DATA IN ONE STRUCTURE FOR USE IN SUBROUTINES
data.KE = KE; data.iK = iK; data.jK = jK; 
data.nelx = nelx; data.nely = nely; data.F = F;
data.freedofs = freedofs; data.edofMat = edofMat;
data.penal = penal; 
data.Emin = Emin; data.E0 = E0; data.volfrac = volfrac;
data.HsH  = HsH;
%% VOLUME CONSTRAINT: g= Ax - volfrac
A = (data.HsH'*repmat(1/(nelx*nely),nelx*nely,1))';
%% LOWER AND UPPER DESIGN BOUNDS
xmin = zeros(nelx*nely,1);
xmax = ones (nelx*nely,1);
%% INITIAL DESIGN
x0 = volfrac(ones(nelx*nely,1));

%% CALL THE OPTIMIZATION SOFTWARE
switch optimizer
    case 'fmincon'
        %% OPTIMIZATION OPTIONS (help fmincon)
        options = optimoptions('fmincon','Algorithm', 'interior-point', ...
            'SpecifyObjectiveGradient',true, ...
            'CheckGradients', false, ...
            'PlotFcn', {@(x,ov,s) PlotX(x,ov,s,data),@optimplotfval}, ...
            'Display','iter', ...
            'HessianApproximation',{'lbfgs',6});
       
        [x,fval,flag,output] = fmincon( ...
            ... % objective function/gradient:
            @(x) fdf(x,data), ...
            ... % initial guess:
            x0,...
            ... % linear inequality constraints: A*x <= volfrac
            A, volfrac, ...
            ... % linear equality constraints: none
            [],[], ...
            ... % lower/upper bounds
            xmin, xmax, ...
            ... % non-linear constraints: none
            [], ...
            ... % finally, optimization options
            options);
        h=get(gcf,'Children');
        g=get(h(5),'children');
        iter_fmincon=g.XData;
        fval_fmincon=g.YData;
        figure(2);clf;
        plot(iter_fmincon,fval_fmincon);
        
    case 'ipopt'
        funcs.objective = @(x) FE(x,data);
        funcs.gradient = @(x) gradient(x,data);
        funcs.constraints = @(x) A*x/volfrac-1;
        funcs.jacobian = @(x) sparse(A/volfrac);
        funcs.jacobianstructure =@() sparse(A);
        
        options.lb=zeros(size(x0));
        options.ub=ones(size(x0));
        options.cu=zeros(size(A,1));
        options.cl=-inf*ones(size(A,1));
        options.ipopt.hessian_approximation = 'limited-memory';
        options.ipopt.limited_memory_max_history = 6;
        % https://coin-or.github.io/Ipopt/OPTIONS.html
        % options.ipopt.derivative_test = 'first-order';
        % options.ipopt.derivative_test_perturbation= 1e-3;
        [x, info] = ipopt(x0,funcs,options);
        [stop]=PlotX(x,[],[],data);
        
    case 'mma'
       error('MMA not implemented yet');
    case 'oc'
        %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        x=x0;
        for loop = 1:300
            [c,dc]=FE(x,data);
            [stop]=PlotX(x,[],[],data);
            dv=A'; fval(loop)=c;
            l1 = 0; l2 = 1e9; move = 0.2;
            while (l2-l1)/(l1+l2) > 1e-3
                lmid = 0.5*(l2+l1);
                xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
                % Investigate feasibility
                if A*xnew-volfrac > 0, l1 = lmid; else l2 = lmid; end
            end
            change = max(abs(xnew(:)-x(:)));
            x = xnew;
            %% PRINT RESULTS
            fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
                A*x,change);
            if change < 0.01
                disp('OC terminated due to max change below 0.01')
                break
            end
            figure(2);hold on
            plot(1:loop,fval,'-r');
        end
    otherwise 
        error('Optimizer not implemented');
end
%------------------------------------------------------------------
%% Objective function and its gradient
function [c,dc] = fdf(x,data)
if(nargout>1)
  [c,dc]=FE(x,data);
else
  [c]=FE(x,data);
end
%------------------------------------------------------------------
function dc = gradient(x,data)
% Plot the design every time the gradient is evaluated
% Multiple function evaluations per iteration is possible for some
% optimizers
PlotX(x,0,0,data); 
% 
[~,dc]=FE(x,data);

%% Plot function
function [stop]=PlotX(x,~,~,data)
%% Filtering
xPhys = data.HsH * x; 
%% PLOT DENSITIES
figure(1);%clf;
colormap(gray); imagesc(reshape(1-xPhys,data.nely,data.nelx)); 
caxis([0 1]); axis ij; axis equal; axis off; drawnow;
stop=false;
%------------------------------------------------------------------
%% perform FE-analysis
function [c,dc]=FE(x,d)
% store some variables between the calls
persistent U x_old xPhys L s ce
if length(x_old) ~= length(x) % only allocate first time
  % pre-allocate memory
  x_old = zeros(size(x));
  U     = zeros(2*(d.nely+1)*(d.nelx+1),1); 
end
% check if design changed since last evaluation
if any(x_old~=x) 
  % need to re-assemble and re-factorize
  x_old = x;
  %% Filtering
  xPhys = d.HsH * x; 
  %% FE-ANALYSIS
  sK = reshape(d.KE(:)*(d.Emin+xPhys(:)'.^d.penal*(d.E0-d.Emin)),...
               64*d.nelx*d.nely,1);
  K = sparse(d.iK,d.jK,sK); K = (K+K')/2;
  % Cholesky factorization
  [L,~,s]=chol(K(d.freedofs,d.freedofs),'lower','vector');
  % Forward/backward substitution
  U(d.freedofs(s))=L'\(L\d.F(d.freedofs(s)));
  %
  ce =  sum((U(d.edofMat)*d.KE).*U(d.edofMat),2);
else
%     disp('Solution not changed');
end
%
% compute outputs
%
%FIXIT: compute objective function (compliance)
 c  =  
if nargout > 1
  %FIXIT: compute sensitivities of compliance
  dc = 

  %% MODIFICATION OF SENSITIVITIES
  dc = d.HsH' * dc;
end
