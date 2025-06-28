function top88_dynamic(nelx,nely,volfrac,penal,rmin,ft,omega_fraction)

% Verify density filter is used (as per hints)
if ft ~= 2
    warning('Hints recommend using density filter (ft=2). Current ft=%d', ft);
    fprintf('Consider running with: top88_dynamic(%d,%d,%.1f,%.1f,%.1f,2,%.1f)\n', ...
        nelx,nely,volfrac,penal,rmin,omega_fraction);
end
%% MATERIAL PROPERTIES
E0 = 100;           % Young's modulus (changed from 1 to 100 as suggested)
Emin = 1e-4;        % Minimum stiffness (as suggested in hints)
nu = 0.3;           % Poisson's ratio
rho0 = 1;           % Unit mass density
rhomin = 1e-8;      % Minimum mass density

%% PREPARE FINITE ELEMENT ANALYSIS
% Stiffness matrix (same as original)
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

% Mass matrix (4-node element, from Jensen paper appendix)
ME = [4/9 0 2/9 0 1/9 0 2/9 0
      0 4/9 0 2/9 0 1/9 0 2/9
      2/9 0 4/9 0 2/9 0 1/9 0
      0 2/9 0 4/9 0 2/9 0 1/9
      1/9 0 2/9 0 4/9 0 2/9 0
      0 1/9 0 2/9 0 4/9 0 2/9
      2/9 0 1/9 0 2/9 0 4/9 0
      0 2/9 0 1/9 0 2/9 0 4/9]/4;

% Element connectivity
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

%% PHASE 1: STATIC OPTIMIZATION TO GET INITIAL DESIGN
fprintf('=== PHASE 1: Static Optimization ===\n');
x = repmat(volfrac,nely,nelx);
xPhys = x;
loop = 0;
change = 1;

% MMA INITIALIZATION
m = 1;
n = nelx*nely;
xmin = zeros(n,1);
xmax = ones(n,1);
xold1 = x(:);
xold2 = x(:);
low = [];
upp = [];
a0 = 1;
a = zeros(m,1);
c = 1000*ones(m,1);
d = zeros(m,1);
% Move limits for MMA (as suggested in hints)
move = 0.2;

% Static optimization loop
while change > 0.01 && loop < 100
  loop = loop + 1;
  
  % FE-ANALYSIS (Static)
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  
  % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS (Static compliance)
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  dv = ones(nely,nelx);
  
  % FILTERING
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  
  % MMA UPDATE
  xval = x(:);
  f0val = c;
  df0dx = dc(:);
  fval = [sum(xPhys(:))/(volfrac*nelx*nely) - 1];
  dfdx = dv(:)'/(volfrac*nelx*nely);
  
  % Apply move limits
  xmin_iter = max(0, x(:) - move);
  xmax_iter = min(1, x(:) + move);
  
  [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
    mmasub(m,n,loop,xval,xmin_iter,xmax_iter,xold1,xold2, ...
    f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
  
  xnew = reshape(xmma,nely,nelx);
  
  if ft == 1
    xPhys = xnew;
  elseif ft == 2
    xPhys(:) = (H*xnew(:))./Hs;
  end
  
  change = max(abs(xnew(:)-x(:)));
  xold2 = xold1;
  xold1 = x(:);
  x = xnew;
  
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
end

%% PHASE 2: FREQUENCY ANALYSIS OF STATIC DESIGN
fprintf('\n=== PHASE 2: Frequency Analysis ===\n');

% Store static design for comparison
x_static = xPhys;

% Assemble mass matrix for static design
sM = reshape(ME(:)*(rhomin+xPhys(:)'.*(rho0-rhomin)),64*nelx*nely,1);
M = sparse(iK,jK,sM); M = (M+M')/2;

% Find first few eigenfrequencies
try
    [phi, lambda] = eigs(K(freedofs,freedofs), M(freedofs,freedofs), 10, 'smallestabs');
    eigenfreqs = sqrt(real(diag(lambda)))/(2*pi);
    fprintf('First 5 eigenfrequencies: %.3f, %.3f, %.3f, %.3f, %.3f Hz\n', eigenfreqs(1:5));
    
    omega1 = sqrt(real(lambda(1,1)));  % First natural frequency (rad/s)
    omega_target = omega_fraction * omega1;  % Target frequency
    
    fprintf('First natural frequency: %.3f rad/s (%.3f Hz)\n', omega1, omega1/(2*pi));
    fprintf('Target frequency (%.1f%% of first): %.3f rad/s (%.3f Hz)\n', ...
        omega_fraction*100, omega_target, omega_target/(2*pi));
    fprintf('*** To change target frequency, modify omega_fraction parameter ***\n');
    fprintf('*** Current: %.1f (try 0.5, 0.8, 0.9, 1.1, 1.5) ***\n', omega_fraction);
catch
    warning('Eigenvalue analysis failed. Using default frequency.');
    omega_target = 10;  % Default target frequency
end

% Plot static design
figure('Position', [100, 100, 1200, 400]);
subplot(1,2,1);
colormap(gray); imagesc(1-x_static); caxis([0 1]); axis equal; axis off;
title('Static Compliance Optimized Design');

% Leave space for dynamic design (will be filled later)
subplot(1,2,2);
text(0.5, 0.5, 'Dynamic optimization in progress...', 'HorizontalAlignment', 'center');
axis off;

% Plot frequency response of static design
omega_range = linspace(0.1*omega1, 1.5*omega1, 100);
response_static = zeros(size(omega_range));

fprintf('Computing frequency response of static design...\n');
for i = 1:length(omega_range)
    omega = omega_range(i);
    S = K - omega^2 * M;
    try
        u_temp = S(freedofs,freedofs) \ F(freedofs);
        response_static(i) = abs(u_temp' * F(freedofs));
    catch
        response_static(i) = inf;  % Near resonance
    end
end

%% PHASE 3: DYNAMIC OPTIMIZATION
fprintf('\n=== PHASE 3: Dynamic Optimization ===\n');

% Reset MMA parameters for dynamic optimization
loop = 0;
change = 1;
xold1 = x(:);
xold2 = x(:);
low = [];
upp = [];
move = 0.2;  % Move limits for dynamic optimization

% Dynamic optimization loop
while change > 0.01 && loop < 100
  loop = loop + 1;
  
  % ASSEMBLE DYNAMIC SYSTEM
  % Stiffness matrix
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  
  % Mass matrix (linear interpolation)
  sM = reshape(ME(:)*(rhomin+xPhys(:)'.*(rho0-rhomin)),64*nelx*nely,1);
  M = sparse(iK,jK,sM); M = (M+M')/2;
  
  % Dynamic stiffness matrix
  S = K - omega_target^2 * M;
  
  % Solve dynamic system
  u_full = zeros(2*(nely+1)*(nelx+1),1);
  u_full(freedofs) = S(freedofs,freedofs) \ F(freedofs);
  
  % DYNAMIC COMPLIANCE OBJECTIVE (Jensen Section 3.2)
  c = abs(u_full' * F);  % |u^T f|
  
  % SENSITIVITY ANALYSIS (Jensen Section 3.2.1 - undamped case)
  alpha = sign(real(u_full' * F));
  lambda_full = -alpha * u_full;
  
  % Element-wise sensitivities
  dc = zeros(nely,nelx);
  for elx = 1:nelx
    for ely = 1:nely
      e = (elx-1)*nely + ely;
      edof = edofMat(e,:);
      
      % Stiffness sensitivity
      dKe = penal * xPhys(ely,elx)^(penal-1) * (E0-Emin);
      
      % Mass sensitivity (linear interpolation)
      dMe = (rho0-rhomin);
      
      % Dynamic stiffness sensitivity: dS/dx = dK/dx - omega^2 * dM/dx
      dSe = dKe * KE - omega_target^2 * dMe * ME;
      
      % Total sensitivity (Jensen Eq. 26, undamped case Eq. 41)
      dc(ely,elx) = -alpha * real(u_full(edof)' * dSe * u_full(edof));
    end
  end
  
  % Volume sensitivity
  dv = ones(nely,nelx);
  
  % FILTERING
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  
  % MMA UPDATE
  xval = x(:);
  f0val = c;
  df0dx = dc(:);
  fval = [sum(xPhys(:))/(volfrac*nelx*nely) - 1];
  dfdx = dv(:)'/(volfrac*nelx*nely);
  
  % Apply move limits
  xmin_iter = max(0, x(:) - move);
  xmax_iter = min(1, x(:) + move);
  
  [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
    mmasub(m,n,loop,xval,xmin_iter,xmax_iter,xold1,xold2, ...
    f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
  
  xnew = reshape(xmma,nely,nelx);
  
  if ft == 1
    xPhys = xnew;
  elseif ft == 2
    xPhys(:) = (H*xnew(:))./Hs;
  end
  
  change = max(abs(xnew(:)-x(:)));
  xold2 = xold1;
  xold1 = x(:);
  x = xnew;
  
  fprintf(' It.:%5i Dyn.Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
  
  % Plot current design every 10 iterations
  if mod(loop,10) == 0 || loop == 1
    % Update the dynamic design plot
    figure(1);
    subplot(1,2,2);
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; 
    title(sprintf('Dynamic Optimized Design (Iter %d)', loop));
    drawnow;
  end
end

%% PHASE 4: FREQUENCY RESPONSE COMPARISON
fprintf('\n=== PHASE 4: Frequency Response Comparison ===\n');

% Compute frequency response of optimized design
response_optimized = zeros(size(omega_range));

% Final system matrices
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K_final = sparse(iK,jK,sK); K_final = (K_final+K_final')/2;
sM = reshape(ME(:)*(rhomin+xPhys(:)'.*(rho0-rhomin)),64*nelx*nely,1);
M_final = sparse(iK,jK,sM); M_final = (M_final+M_final')/2;

fprintf('Computing frequency response of optimized design...\n');
for i = 1:length(omega_range)
    omega = omega_range(i);
    S = K_final - omega^2 * M_final;
    try
        u_temp = S(freedofs,freedofs) \ F(freedofs);
        response_optimized(i) = abs(u_temp' * F(freedofs));
    catch
        response_optimized(i) = inf;
    end
end

% Plot comparison
figure;
semilogy(omega_range/(2*pi), response_static, 'b-', 'LineWidth', 2); hold on;
semilogy(omega_range/(2*pi), response_optimized, 'r-', 'LineWidth', 2);
xline(omega_target/(2*pi), 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
xlabel('Frequency [Hz]');
ylabel('Dynamic Compliance |u^T f|');
legend('Static Design', 'Dynamic Optimized', 'Target Frequency', 'Location', 'best');
title('Frequency Response Comparison');
grid on;

fprintf('\n=== OPTIMIZATION COMPLETED ===\n');
fprintf('Target frequency: %.3f Hz (%.1f%% of first resonance)\n', ...
    omega_target/(2*pi), omega_fraction*100);
fprintf('Static design response at target: %.3e\n', ...
    interp1(omega_range, response_static, omega_target));
fprintf('Optimized design response at target: %.3e\n', ...
    interp1(omega_range, response_optimized, omega_target));

% Calculate improvement
improvement = interp1(omega_range, response_static, omega_target) / ...
              interp1(omega_range, response_optimized, omega_target);
fprintf('Improvement factor: %.2fx\n', improvement);

fprintf('\n*** TO CHANGE TARGET FREQUENCY ***\n');
fprintf('Modify the omega_fraction parameter:\n');
fprintf('  omega_fraction = 0.5  → 50%% of ω₁ (well below resonance)\n');
fprintf('  omega_fraction = 0.8  → 80%% of ω₁ (approaching resonance)\n');
fprintf('  omega_fraction = 0.9  → 90%% of ω₁ (close to resonance) [current]\n');
fprintf('  omega_fraction = 1.0  → 100%% of ω₁ (at resonance - challenging!)\n');
fprintf('  omega_fraction = 1.1  → 110%% of ω₁ (above resonance)\n');
fprintf('  omega_fraction = 1.5  → 150%% of ω₁ (well above resonance)\n');

%% FINAL PLOT
figure;
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off;
title('Final Dynamic Topology Optimized Design');

end