%% PBM‑8 – Time‑harmonic MBB beam workflow 
%   Jensen (2017) – Problem 8 helper script
%   ------------------------------------------------------------
%   1)  Run a *static* compliance optimisation (ω = 0) and keep the
%       converged density field `x_static`.
%   2)  Evaluate the frequency response of this design over a user‑defined
%       band of angular frequencies and plot the objective ‖uᵀf‖ on a log
%       y‑scale.
%   3)  Extract the *first* resonance frequency ω₁ (the first peak in the
%       response) and define an optimisation target ωₜ = 0.9 ω₁.
%   4)  Re‑optimise the beam for *dynamic* compliance at ωₜ, starting from
%       `x_static`.
%   5)  Compare the frequency responses of the initial and the optimised
%       design in a single plot.
%   ------------------------------------------------------------
%   E. Andreassen et al. (2010) 88‑line code + small extensions by you :-)
%   ------------------------------------------------------------
%% House‑keeping
clear;  clc;  close all;

%% ---------- USER SETTINGS ------------------------------------
nelx   = 100;            % number of elements in x–direction
nely   =  50;            % number of elements in y–direction
volfrac = 0.5;           % prescribed volume fraction
penal   = 3.0;           % SIMP penalty
rmin    = 1.2;           % filter radius
ft      = 2;             % 1=density, 2=stiffness filter
method  = "mma";         % {"oc","mma"} optimisation update

% Frequency sweep (rad s⁻¹).  Avoid zero exactly → start at 1e‑6 :
omega_sweep = logspace(-1, 1, 100);   % ~10⁻⁶ … 10³   (adjust if needed)   % ~10⁻⁶ … 10³   (adjust if needed)
%% -------------------------------------------------------------

%% 1) Static compliance optimisation (ω = 0)
[x_static, ~] = top88_pbm8(nelx,nely,volfrac,penal,rmin,ft,method,0,[]);

%% 2) Frequency response of the static design
c_static = arrayfun(@(om) dynCompliance(x_static,om, ...
    nelx,nely,volfrac,penal,rmin,ft), omega_sweep);

figure(1);
semilogy(omega_sweep, c_static,'b','LineWidth',1.2); hold on;
xlabel('\omega  [rad/s]'); ylabel('|u^{T}f|  (dynamic compliance)');
set(gca,'XScale','log'); grid on;

%% 3) Locate first resonance (first *local* maximum)
[~,idx] = max(c_static);              % crude – first dominant peak
omega1  = omega_sweep(idx);
omega_t = 0.9*omega1;
fprintf('\nFirst resonance ω₁ ≈ %.4g   →  optimisation target ω_t = %.4g\n', omega1, omega_t);

%% 4) Dynamic compliance optimisation at ω_t (warm‑start with x_static)
[x_dyn, ~] = top88_pbm8(nelx,nely,volfrac,penal,rmin,ft,method,omega_t,x_static);

%% 5) Frequency response of the optimised design
c_dyn = arrayfun(@(om) dynCompliance(x_dyn,om, ...
    nelx,nely,volfrac,penal,rmin,ft), omega_sweep);

semilogy(omega_sweep, c_dyn,'r','LineWidth',1.2);
legend({'initial (ω = 0)','optimised (ω_t)'},'Location','best');

title(sprintf('Frequency response – initial vs optimised (ω_t = %.3g rad/s)',omega_t));

%% -------- OPTIONAL: experiment with other target frequencies --------
%  Uncomment and set a different ω_t (e.g. ω_t = 1.05*ω₁) to inspect
%  behaviour above resonance.
% omega_t = 1.05*omega1;  % example
% [x_dyn2,~] = top88_pbm8(nelx,nely,volfrac,penal,rmin,ft,method,omega_t,x_static);
% c_dyn2 = arrayfun(@(om) dynCompliance(x_dyn2,om, ...
%        nelx,nely,volfrac,penal,rmin,ft), omega_sweep);
% semilogy(omega_sweep, c_dyn2,'g','LineWidth',1.2);
% legend({'initial','optimised 0.9ω₁','optimised 1.05ω₁'});

%% ---------------- HELPER FUNCTION -----------------------------------
function c = dynCompliance(xPhys,omega,nelx,nely,volfrac,penal,rmin,ft)
%   One‑off dynamic compliance evaluation for a *frozen* density field.
%   Replicates the FE analysis block of top88_pbm8 without design updates.
%   (Only the minimum material density rhomin is hard‑coded.)

%% --- material & FE constants (duplicated from top88_pbm8) ---
E0 = 1; Emin = 1e-4; nu = 0.3; rhomin = 1e-8;
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE  = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
% lumped 8×8 element mass matrix (identical to code in top88_pbm8)
m0  = [4/9 0 2/9 0 1/9 0 2/9 0; ...
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

% assemble K & M for the *given* density field
sK = reshape(KE(:)*(Emin + xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K  = sparse(kron(edofMat,ones(8,1))',kron(edofMat,ones(1,8))',sK);
K  = (K+K')./2;

sM = reshape(m0(:)*(rhomin + (1-rhomin)*xPhys(:)'),64*nelx*nely,1);
M  = sparse(kron(edofMat,ones(8,1))',kron(edofMat,ones(1,8))',sM);
M  = (M+M')./2;

S  = K - omega^2*M;

% boundary conditions (half MBB – identical to top88_pbm8)
F = sparse(2,1,- 1,2*(nely+1)*(nelx+1),1);
fixeddofs = union(1:2:2*(nely+1), 2*(nelx+1)*(nely+1));
freedofs  = setdiff(1:2*(nely+1)*(nelx+1),fixeddofs);

U = zeros(2*(nely+1)*(nelx+1),1);
U(freedofs) = S(freedofs,freedofs) \ F(freedofs);

% dynamic compliance (objective) – same formula as in optimiser
ceKU = sum((U(edofMat)*KE).*U(edofMat),2);
ceMU = sum((U(edofMat)*m0).*U(edofMat),2);
ceKU = reshape(ceKU,nely,nelx);
ceMU = reshape(ceMU,nely,nelx);

c = norm( sum(sum((Emin + xPhys.^penal*(E0-Emin)).*ceKU)) - omega^2 * sum(sum((rhomin + (1-rhomin)*xPhys).*ceMU)) );
end
