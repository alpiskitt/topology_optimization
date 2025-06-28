function top88_dynamic_fd_check(nelx,nely,volfrac,penal,rmin,ft,omega_fraction)

% Verify density filter is used (as per hints)
if ft ~= 2
    warning('Hints recommend using density filter (ft=2). Current ft=%d', ft);
    fprintf('Consider running with: top88_dynamic_fd_check(%d,%d,%.1f,%.1f,%.1f,2,%.1f)\n', ...
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

%% INITIALIZE WITH UNIFORM DESIGN
x = repmat(volfrac,nely,nelx);
if ft == 1
    xPhys = x;
elseif ft == 2
    xPhys(:) = (H*x(:))./Hs;
end

% Set target frequency
omega_target = omega_fraction * 10;  % Use a reasonable default frequency

%% FINITE DIFFERENCE CHECK
fprintf('\n=== FINITE DIFFERENCE CHECK ===\n');

% Function to compute dynamic compliance
compute_dynamic_compliance = @(x_input) compute_compliance(x_input, nelx, nely, KE, ME, ...
    E0, Emin, rho0, rhomin, penal, omega_target, edofMat, iK, jK, freedofs, F, H, Hs, ft);

% Current design point
x0 = x;
[c0, dc_analytical] = compute_dynamic_compliance(x0);

% Finite difference parameters
h_values = 10.^(-1:-1:-6);  % Step sizes from 1e-1 to 1e-12
n_elements = min(10, nelx*nely);  % Check first 10 elements (or all if fewer)
errors = zeros(length(h_values), n_elements);

fprintf('Computing finite differences...\n');
fprintf('Analytical objective: %.6e\n', c0);

% Perform finite difference check for selected elements
for i = 1:n_elements
    % Convert linear index to subscripts
    [ely, elx] = ind2sub([nely, nelx], i);
    
    fprintf('Element (%d,%d): dc_analytical = %.6e\n', elx, ely, dc_analytical(ely,elx));
    
    for j = 1:length(h_values)
        h = h_values(j);
        
        % Forward difference
        x_plus = x0;
        x_plus(ely,elx) = x_plus(ely,elx) + h;
        c_plus = compute_dynamic_compliance(x_plus);
        
        % Backward difference
        x_minus = x0;
        x_minus(ely,elx) = x_minus(ely,elx) - h;
        c_minus = compute_dynamic_compliance(x_minus);
        
        % Central difference
        dc_numerical = (c_plus - c_minus) / (2 * h);
        
        % Compute relative error
        if abs(dc_analytical(ely,elx)) > 1e-12
            errors(j, i) = abs(dc_numerical - dc_analytical(ely,elx)) / abs(dc_analytical(ely,elx));
        else
            errors(j, i) = abs(dc_numerical - dc_analytical(ely,elx));
        end
        
        if j <= 3  % Print first few for debugging
            fprintf('  h=%.1e: dc_numerical=%.6e, error=%.2e\n', h, dc_numerical, errors(j,i));
        end
    end
end

%% PLOT RESULTS
fprintf('\nPlotting finite difference check results...\n');

figure('Position', [100, 100, 1200, 800]);

% Plot 1: Error vs step size for all checked elements
subplot(2,2,1);
loglog(h_values, errors, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
% Add reference lines
loglog(h_values, h_values.^2, 'k--', 'LineWidth', 2, 'DisplayName', 'h^2 (ideal)');
loglog(h_values, ones(size(h_values))*1e-15, 'r--', 'LineWidth', 2, 'DisplayName', 'Machine precision');
xlabel('Step size h');
ylabel('Relative error');
title('Finite Difference Check: All Elements');
legend('Location', 'best');
grid on;
xlim([min(h_values), max(h_values)]);
ylim([1e-16, 1e2]);

% Plot 2: Best errors for each element
subplot(2,2,2);
min_errors = min(errors, [], 1);
optimal_h = h_values(1) * ones(size(min_errors));
for i = 1:n_elements
    [~, idx] = min(errors(:,i));
    optimal_h(i) = h_values(idx);
end

semilogy(1:n_elements, min_errors, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Element number');
ylabel('Minimum relative error');
title('Best Achievable Error per Element');
grid on;

% Plot 3: Optimal step size for each element
subplot(2,2,3);
loglog(1:n_elements, optimal_h, 'bs-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Element number');
ylabel('Optimal step size');
title('Optimal Step Size per Element');
grid on;

% Plot 4: Error distribution
subplot(2,2,4);
histogram(log10(min_errors), 'BinWidth', 0.5, 'FaceColor', 'blue', 'EdgeColor', 'black');
xlabel('log_{10}(minimum relative error)');
ylabel('Count');
title('Distribution of Minimum Errors');
grid on;

% Summary statistics
fprintf('\n=== FINITE DIFFERENCE CHECK SUMMARY ===\n');
fprintf('Number of elements checked: %d\n', n_elements);
fprintf('Best overall error: %.2e\n', min(min_errors));
fprintf('Worst error: %.2e\n', max(min_errors));
fprintf('Median error: %.2e\n', median(min_errors));
fprintf('Mean error: %.2e\n', mean(min_errors));

% Check if gradients are correct
if max(min_errors) < 1e-6
    fprintf('✓ GRADIENTS APPEAR CORRECT (all errors < 1e-6)\n');
elseif max(min_errors) < 1e-3
    fprintf('⚠ GRADIENTS MOSTLY CORRECT (all errors < 1e-3)\n');
else
    fprintf('✗ GRADIENT ERRORS DETECTED (some errors > 1e-3)\n');
end

% Plot design
figure('Position', [100, 950, 600, 300]);
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off;
title('Design Used for Finite Difference Check');

end

%% HELPER FUNCTION TO COMPUTE DYNAMIC COMPLIANCE AND SENSITIVITIES
function [c, dc] = compute_compliance(x, nelx, nely, KE, ME, E0, Emin, rho0, rhomin, ...
    penal, omega_target, edofMat, iK, jK, freedofs, F, H, Hs, ft)

% Apply filter
if ft == 1
    xPhys = x;
elseif ft == 2
    xPhys = reshape((H*x(:))./Hs, nely, nelx);
end

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
try
    u_full(freedofs) = S(freedofs,freedofs) \ F(freedofs);
catch
    warning('Singular matrix encountered, using pseudoinverse');
    u_full(freedofs) = pinv(full(S(freedofs,freedofs))) * F(freedofs);
end

% DYNAMIC COMPLIANCE OBJECTIVE
c = abs(u_full' * F);

% Only compute sensitivities if requested
if nargout > 1
    % SENSITIVITY ANALYSIS (Jensen Section 3.2.1 - undamped case)
    alpha = sign(real(u_full' * F));
    
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
    
    % Apply filter to sensitivities
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
    end
end

end