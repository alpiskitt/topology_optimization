function topS140_load(nelx,nely,volfrac,penal,rmin,ft,ftBC,eta,beta,move,pnorm,maxit,type)
% ---------------------------- PRE. 1) MATERIAL AND CONTINUATION PARAMETERS
E0 = 1;                                                                     % Young modulus of solid
Emin = 1e-9;                                                                % Young modulus of "void"
nu = 0.3;                                                                   % Poisson ratio
penalCnt = { 3,  5, 50, 0 };                                                % continuation penal
betaCnt  = { 2,  16, 75, 1 };                                               % continuation beta
if ftBC == 'N', bcF = 'symmetric'; else, bcF = 0; end                       % filter BC selector
% ----------------------------------------- PRE. 2) DISCRETIZATION FEATURES
nEl = nelx * nely;                                                          % number of elements
nodeNrs = int32( reshape( 1 : (1 + nelx) * (1 + nely), 1+nely, 1+nelx ) );  % nodes numbers
cVec = reshape( 2 * nodeNrs( 1 : end - 1, 1 : end - 1 ) + 1, nEl, 1 );
cMat = cVec + int32( [ 0, 1, 2 * nely + [ 2, 3, 0, 1 ], -2, -1 ] );         % connectivity matrix
nDof = ( 1 + nely ) * ( 1 + nelx ) * 2;                                     % total number of DOFs
[ sI, sII ] = deal( [ ] );
for j = 1 : 8
    sI = cat( 2, sI, j : 8 );
    sII = cat( 2, sII, repmat( j, 1, 8 - j + 1 ) );
end
[ iK , jK ] = deal( cMat( :,  sI )', cMat( :, sII )' );
Iar = sort( [ iK( : ), jK( : ) ], 2, 'descend' ); clear iK jK
c1 = [12;3;-6;-3;-6;-3;0;3;12;3;0;-3;-6;-3;-6;12;-3;0;-3;-6;3;12;3;...
    -6;3;-6;12;3;-6;-3;12;3;0;12;-3;12];
c2 = [-4;3;-2;9;2;-3;4;-9;-4;-9;4;-3;2;9;-2;-4;-3;4;9;2;3;-4;-9;-2;...
    3;2;-4;3;-2;9;-4;-9;4;-4;-3;-4];
Ke = 1/(1-nu^2)/24*( c1 + nu .* c2 );                                       % lower sym. part
Ke0( tril( ones( 8 ) ) == 1 ) = Ke';
Ke0 = reshape( Ke0, 8, 8 );
Ke0 = Ke0 + Ke0' - diag( diag( Ke0 ) );                                     % recover full matrix
% ------------------------------------ PRE. 3) SUPPORTS AND PASSIVE DOMAINS
nrlo = union(nodeNrs(round(nely/2):end,1),nodeNrs(round(nely/2):end,end));
fixed = union(2*nrlo,2*nrlo-1);                                             % fix in both directions
[ pasS, pasV ] = deal(1:nely:nelx*nely,[]);
free = setdiff( 1 : nDof, fixed );                                          % set of free DOFs
act = setdiff( (1 : nEl )', union( pasS, pasV ) );                          % set of active d.v.
% ----------------------------------------------------------- PRE. 4) LOADS
lcDof = 2 * nodeNrs(1,:);
% ------------------------------------------------ PRE. 5) STOCHASTIC MODEL
rng('default');                                                             % reset random parameter
com0 = 100;                                                                 % initial guess (scaling)
y_all = load('Data.mat','r'); y_all = round(nelx*y_all.r)+1; 
y = unique(y_all); n_disc = size(y,2);                                      % integration points
dist_w = sum(y_all' == y)/size(y_all,2);                                    % weights quadrature rule
switch type
    case 'distribution' % sample according to distribution
        X = y_all(randi(numel(y_all),1,maxit));
    case 'uniform'  % sample uniformly
        X = y(randi(n_disc,1,maxit));
end
maxsmpl = 2000;                                                             % max stored samples
[x_birth,x_ind,leavers] = deal(1:maxsmpl,1:maxsmpl,1);
y_weight = 5*volfrac*sqrt(nEl);  % weighting norm
%keyboard
y_diff = pdist2(y',X');
y_diff = y_diff/max(max(y_diff,1e-10),[],'all');                            % normalization
X = int32(X);
[Gra,DesH,ComH] = deal( zeros(nEl,maxsmpl), zeros(nEl,maxsmpl), zeros(1,maxsmpl) );
% --------------------------------------- PRE. 6) DEFINE IMPLICIT FUNCTIONS
prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v(:)-eta)))./...
    (tanh(beta*eta)+tanh(beta*(1-eta)));                                    % projection
deta = @(v,eta,beta) - beta * csch( beta ) .* sech( beta * ( v( : ) - eta ) ).^2 .* ...
    sinh( v( : ) * beta ) .* sinh( ( 1 - v( : ) ) * beta );
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));
cnt = @(v,vCnt,l) v+(l>=vCnt{1})*(v<vCnt{2})*(mod(l,vCnt{3})==0)*vCnt{4};   % apply continuation
% -------------------------------------------------- PRE. 7) PREPARE FILTER
[dy,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
h = max( 0, rmin - sqrt( dx.^2 + dy.^2 ) );                                 % conv. kernel
Hs = imfilter( ones( nely, nelx ), h, bcF );                                % filter weights
dHs = Hs;
% ------------------------ PRE. 8) ALLOCATE AND INITIALIZE OTHER PARAMETERS
[ x, dsK, dV ] = deal( zeros( nEl, 1 ) );                                   % initialize vectors
dV( act, 1 ) = 1/nEl/volfrac;                                               % derivative of volume
x( act ) = ( volfrac*( nEl - length(pasV) ) - length(pasS) )/length( act ); % volume fraction
x( pasS ) = 1;                                                              % set x = 1 on pasS set
[ xPhys, loop, U ] = deal( x, 0, zeros( nDof, 1 ) );                        % it. counter, U
% ================================================= START OPTIMIZATION LOOP
while loop < maxit 
    loop = loop + 1;                                                        % update iter. counter
    % --------- RL. 1) COMPUTE PHYSICAL DENSITY FIELD (AND ETA IF PROJECT.)
    xTilde = imfilter( reshape( x, nely, nelx ), h, bcF ) ./ Hs;
    xPhys( act ) = xTilde( act );  
    if ft > 1                              % compute optimal eta* with Newton
        f = ( mean( prj( xPhys, eta, beta ) ) - volfrac ) * ( ft == 3 );    % function (volume)
        while abs( f ) > 1e-6           % Newton process for finding opt. eta
            eta = eta - f / mean( deta( xPhys( : ), eta, beta ) );
            f = mean( prj( xPhys, eta, beta ) ) - volfrac;
        end
        dHs = Hs ./ reshape( dprj( xTilde, eta, beta ), nely, nelx );       % modification sensitivity
        xPhys = prj( xPhys, eta, beta );                                    % projected (phys.) field
    end
    % ------------------------ RL. 2) SETUP AND SOLVE EQUILIBRIUM EQUATIONS
    F = fsparse( lcDof(X(:,loop)), 1, -1, [ nDof, 1 ] );
    sK = ( Emin + xPhys.^penal * ( E0 - Emin ) );                           % stiffness interpolation
    dsK( act ) = -penal * ( E0 - Emin ) * xPhys( act ) .^ ( penal - 1 );
    sK = reshape( Ke( : ) * sK', length( Ke ) * nEl, 1 );
    K = fsparse( Iar( :, 1 ), Iar( :, 2 ), sK, [ nDof, nDof ] );            % stiffness matrix
    U(free) = decomposition( K( free, free ), 'chol','lower' ) \ F(free);   % solve equilibrium system
    % ---------------------------------------- RL. 3) COMPUTE SENSITIVITIES
    dV0 = imfilter( reshape( dV, nely, nelx ) ./ dHs, h, bcF );
    dc = reshape(dsK .* sum( (U( cMat ) * Ke0 ) .* U( cMat ), 2 ), nely, nelx);
    dc = imfilter( dc ./ dHs, h, bcF );
    % -------------------- RL. 4) SAMPLE MANAGEMENT AND INTEGRATION WEIGHTS
    ulim = min(maxsmpl,loop);
    [Gra(:,leavers), ComH(:,leavers), DesH(:,leavers)] = deal( dc(:), F'*U, xPhys );
    [~,csw] = min( vecnorm(xPhys-DesH(:,1:ulim),2,1) + y_weight*y_diff(:,x_ind(1:ulim)), [], 2);
    weights = sum((csw==1:ulim).*dist_w');                                  % integration weights
    if (loop+1) <= maxsmpl
        leavers = loop+1;
    else
        ind_can = find(weights-min(weights)<1e-8);                          % small weights
        [~,iind] = min(x_birth(ind_can));                                   % oldest of these samples
        leavers = ind_can(iind);
        [x_ind(leavers), x_birth(leavers)] = deal(loop+1,loop+1);
    end
    % ------------------------------ RL. 5) NEAREST NEIGHBOR APPROXIMATIONS
    Compl = sum(weights.*ComH(:,1:ulim),2);Cp = (sum(weights.*(ComH(:,1:ulim).^pnorm),2))^(1/pnorm);
    if mod(loop,25) == 0
        com0 = Compl;                                                       % adjust normalization
    end
    dc = 1/com0*sum(weights.*(((ComH(:,1:ulim)/com0).^(pnorm-1)).*Gra(:,1:ulim)),2);
    % --------------- RL. 6) UPDATE DESIGN VARIABLES AND APPLY CONTINUATION
    xT = x( act );
    [ xU, xL ] = deal( xT + move, xT - move );                              % current bounds
    ocP = xT .* real( sqrt( -dc( act ) ./ dV0( act ) ) );
    l = [ 0, 10*(mean( ocP ) / volfrac) ];                                  % initial estimate for LM
    while ( l( 2 ) - l( 1 ) ) / ( l( 2 ) + l( 1 ) ) > 1e-8  && l(2) > 1e-40 % OC resizing rule
        lmid = 0.5 * ( l( 1 ) + l( 2 ) );
        x( act ) = max( max( min( min( ocP / lmid, xU ), 1 ), xL ), 0 );
        xTilde = imfilter( reshape( x, nely, nelx ), h, bcF ) ./ Hs;
        xPhys( act ) = xTilde( act );
        xPhys = prj( xPhys, eta, beta );
        if mean( xPhys ) > volfrac, l( 1 ) = lmid; else, l( 2 ) = lmid; end
    end
    [penal,beta] = deal(cnt(penal,penalCnt,loop), cnt(beta,betaCnt,loop));  % apply conitnuation
    % ------------------------ RL. 7) PRINT CURRENT RESULTS AND PLOT DESIGN
    fprintf( 'It.:%5i C:%7.4f Cp:%7.4f V:%7.3f penal:%7.2f beta:%7.1f eta:%7.2f \n', ...
        loop, Compl, Cp, mean(xPhys), penal, beta, eta);
    colormap( gray ); imagesc( 1 - reshape( xPhys, nely, nelx ) );
    caxis([0 1]); axis equal off; drawnow
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by A. Uihlein[1], O. Sigmund[2] and M. Stingl[1]  %
% [1] Dept. of Mathematics, Friedrich-Alexander-Universitaet              %
% Erlangen-Nuernberg, 90158 Erlangen                                      %
% [2] Dept. of Solid Mechanics, Technical University of Denmark,          %
% 2800 Lyngby (DK)                                                        %
% Please send your comments to: andrian.uihlein@fau.de                    %
%                                                                         %
% The code is intended for educational purposes and theoretical details   %
% are discussed in the paper Uihlein, A.; Sigmund, O.; Stingl, M. - A 140 %
% line MATLAB code for topology optimization problems with probabilistic  %
% parameters                                                              %
%                                                                         %
% The code as well as a postscript version of the paper can be            %
% downloaded from the web-site: http://www.topopt.dtu.dk                  %
%                                                                         %
% Disclaimer:                                                             %
% The authors reserves all rights but do not guarantee that the code is   %
% free from errors. Furthermore, we shall not be liable in any event      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%