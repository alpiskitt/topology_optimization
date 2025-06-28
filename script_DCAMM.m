
%% Input parameters
nelx = 300; % size X
nely = 100; % size Y
volfrac = 0.3; % volume fraction for constraint
penal = 3; % penalization factor
rmin = 3.2; % filter radius
ft = 2; % filter type
ftBC = 'N'; % filtering boundary conditions 
eta = 0.5; % threshold level
beta = 2; % threshold strength
pnorm = 1; 
move = 1e-2;

%% optimize with nearest neighbor, sample uniformly
close all
topS140_load(nelx,nely,volfrac,penal,rmin,ft,ftBC,eta,beta,move,pnorm,1200,'uniform');
