%% PBM‑8 – Time‑harmonic MBB beam workflow (Jensen 2017, Problem 8)
%   Enhanced version – prettier frequency‑response visualisation
%   -----------------------------------------------------------------------
%   1) Static compliance optimisation (ω = 0) → density field x_static
%   2) Sweep frequencies and plot |uᵀf| for the static design
%   3) Extract first resonance ω₁ and choose ωₜ = 0.9 ω₁
%   4) Dynamic compliance optimisation at ωₜ (warm‑start = x_static)
%   5) Compare the two frequency responses in one semilog‑plot
%   -----------------------------------------------------------------------
%   Author: <you> – built on the 88‑line TOP88 code (Andreassen et al., 2010)
%   -----------------------------------------------------------------------

%% House‑keeping
clear;  clc;  close all;

%% ---------------- USER SETTINGS ---------------------------------------
nelx     = 100;      % elements in x
nely     =  50;      % elements in y
volfrac  = 0.50;     % volume fraction
penal    = 3.0;      % SIMP penalty
rmin     = 1.2;      % filter radius
ft       = 2;        % 1 = density filter, 2 = stiffness filter
method   = "mma";    % {"oc","mma"}

% Frequency sweep (rad s⁻¹) – avoid zero exactly
omega_sweep = logspace(-1, 1, 100);   % 10⁻¹ … 10¹ rad s⁻¹
freq_Hz     = omega_sweep/(2*pi);      % for nicer axis labelling

%% ---------------- 1) STATIC OPTIMISATION -----------------------------
[x_static, ~] = top88_pbm8(nelx,nely,volfrac,penal,rmin,ft,method,0,[]);

%% ---------------- 2) FREQUENCY RESPONSE (static design) --------------
fprintf('\n▷  Evaluating frequency response of STATIC design …\n');
c_static = arrayfun(@(om) dynCompliance(x_static,om, ...
                          nelx,nely,volfrac,penal,rmin,ft), omega_sweep);

%% ---------------- 3) FIRST RESONANCE & TARGET ωₜ ---------------------
[~,idx_pk] = max(c_static);          % first dominant peak ≈ ω₁
omega1  = omega_sweep(idx_pk);
omega_t = 0.9*omega1;
fprintf('   → First resonance  ω₁ ≈ %.4g  rad/s  (%.4g Hz)\n',omega1,omega1/2/pi);
fprintf('   → Target frequency ωₜ = 0.9 ω₁ = %.4g  rad/s\n',omega_t);

%% ---------------- 4) DYNAMIC OPTIMISATION AT ωₜ ----------------------
[x_dyn, ~] = top88_pbm8(nelx,nely,volfrac,penal,rmin,ft,method,omega_t,x_static);

%% ---------------- 5) FREQUENCY RESPONSE (optimised design) -----------
fprintf('\n▷  Evaluating frequency response of OPTIMISED design …\n');
c_dyn = arrayfun(@(om) dynCompliance(x_dyn,om, ...
                      nelx,nely,volfrac,penal,rmin,ft), omega_sweep);

% Optional light smoothing to remove tiny numerical spikes
smoothN     = 3;                            % set to 1 to disable
c_static_pl = movmean(c_static,smoothN,'Endpoints','shrink');
c_dyn_pl    = movmean(c_dyn,   smoothN,'Endpoints','shrink');

%% ----------- PLOT FREQUENCY‑RESPONSE COMPARISON ----------------------
fig = figure('Name','Frequency response comparison','Color','w');

% 1) Static design – blue
semilogy(freq_Hz, c_static_pl, 'Color',[0 0.4470 0.7410], ...
         'LineWidth',2, 'DisplayName','Static design'); hold on;

% 2) Optimised design – red
semilogy(freq_Hz, c_dyn_pl,    'Color',[0.8500 0.3250 0.0980], ...
         'LineWidth',2, 'DisplayName','Optimised design');

% 3) Dashed vertical line at target frequency
xline(omega_t/(2*pi), '--k', 'LineWidth',1.5, 'Alpha',0.7, ...
      'DisplayName','Target frequency');

% Axis cosmetics
set(gca,'XScale','log','YScale','log');
grid on; grid minor;

xlabel('Frequency  f  [Hz]');
ylabel('|u^{T}f|   (dynamic compliance)');

title(sprintf('Frequency response comparison  (ω_t = %.3g  rad/s,  %.3g Hz)', ...
              omega_t, omega_t/2/pi));

legend('Location','best');

% Zoom into resonance region for clarity
xlim([0.8*omega1 1.4*omega1]/(2*pi));
ylim([min([c_static_pl;c_dyn_pl])*0.5 , max([c_static_pl;c_dyn_pl])*2]);

%% ---------------- IMPROVEMENT FACTOR @ ωₜ ----------------------------
C_static_t = interp1(omega_sweep, c_static, omega_t);
C_dyn_t    = interp1(omega_sweep, c_dyn,    omega_t);
improve    = C_static_t / C_dyn_t;

fprintf('\nDynamic compliance at ωₜ:\n');
fprintf('   Static design   : %.3e\n', C_static_t);
fprintf('   Optimised design: %.3e\n', C_dyn_t);
fprintf('   Improvement     : %.2f×\n', improve);

%% ---------------- FINAL DESIGN PLOT ----------------------------------
figure('Name','Final dynamic topology');
colormap(gray); imagesc(1-x_dyn); caxis([0 1]); axis equal off;
colorbar off;
title('Final dynamic‑compliance optimised topology');

%% ---------------- HELPER FUNCTION ------------------------------------
function c = dynCompliance(xPhys,omega,nelx,nely,volfrac,penal,rmin,ft)
%   Single‑shot dynamic compliance for a *fixed* density field xPhys.
%   Replicates the FE analysis of top88_pbm8 without design updates.

%% ---- material & FE constants (duplicated) ----
E0 = 1; Emin = 1e-4; nu = 0.3; rhomin = 1e-8;
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE  = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
% consistent/lumped 8×8 element mass
m0 = [4/9 0 2/9 0 1/9 0 2/9 0; ...
      0 4/9 0 2/9 0 1/9 0 2/9; ...
      2/9 0 4/9 0 2/9 0 1/9 0; ...
      0 2/9 0 4/9 0 2/9 0 1/9; ...
      1/9 0 2/9 0 4/9 0 2/9 0; ...
      0 1/9 0 2/9 0 4/9 0 2/9; ...
      2/9 0 1/9 0 2/9 0 4/9 0; ...
      0 2/9 0 1/9 0 2/9 0 4/9]/4/(nelx*nely);

nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);

iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

%   Assemble stiffness & mass
sK = reshape(KE(:)*(Emin + xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K  = sparse(iK,jK,sK); K = (K+K')./2;
sM = reshape(m0(:)*(rhomin + (1-rhomin)*xPhys(:)'),64*nelx*nely,1);
M  = sparse(iK,jK,sM); M = (M+M')./2;

S  = K - omega^2*M;

%   Loads & BCs (half MBB)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
fixeddofs = union(1:2:2*(nely+1), 2*(nelx+1)*(nely+1));
freedofs  = setdiff(1:2*(nely+1)*(nelx+1), fixeddofs);

U = zeros(2*(nely+1)*(nelx+1),1);
U(freedofs) = S(freedofs,freedofs) \ F(freedofs);

%   Dynamic compliance
ceKU = sum((U(edofMat)*KE).*U(edofMat),2);
ceMU = sum((U(edofMat)*m0).*U(edofMat),2);
ceKU = reshape(ceKU,nely,nelx);
ceMU = reshape(ceMU,nely,nelx);

c = norm( sum(sum((Emin + xPhys.^penal*(E0-Emin)).*ceKU)) ...
        - omega^2*sum(sum((rhomin + (1-rhomin)*xPhys).*ceMU)) );
end
