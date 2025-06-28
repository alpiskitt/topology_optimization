function topS_dynamic(nelx,nely,volfrac,penal,rmin,ft,ftBC,eta,beta,move,pnorm,maxit)
rng default
prjType = 'eta';prjSel = 1;
wInt = [0, .2];
intF    = @( x, p, P0, Pmin )      Pmin + ( P0 - Pmin ) .* x.^p;
dintF   = @( x, p, P0, Pmin ) p * ( P0 - Pmin ) .* x.^( p - 1 );
penalR  = 1;
%% ========================================================================
% MATERIAL AND CONTINUATION PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% ----- store initial values of parameters that are updated by continuation
[ E0, Emin, R0, Rmin, nu ] = deal( 1, 1e-9, 1, 1e-8, 0.3 );                % material properties
[dampK, dampM] = deal(1e-1, 0*1e-1);                                       % Damping coefficients
penalCnt = { 3,  5, 50, 0*.25 };                                           % continuation scheme on penal
betaCnt  = { 2,  16, 50,   1 };                                            % continuation scheme on beta
cnt = @(v,vCn,l) v+(l>=vCn{1}).*(v<vCn{2}).*(mod(l,vCn{3})==0).*vCn{4};    % implicit function for continuation
% DEFINE DISCRETIZATION FEATURES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nEl   = nelx * nely;                                                       % number of elements
[ nodeNrs, nDof, cMat, Ke, Ke0, Me, Me0, Iar] = discQ4( nely, nelx, nu );
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% make complex (Structural / stiffness proportional or mass proportional damping
Ke0 = (1+dampK*1i)*Ke0;
Ke  = (1+dampK*1i)*Ke ;
Me0 = (1+dampM*1i)*Me0;
Me  = (1+dampM*1i)*Me ;
% PREPARE FILTER AND PROJECTION OPERATORS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if ftBC == 'N', bcF = 'symmetric'; else, bcF = 0; end                      % filter BC selector
[dy,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
h = max( 0, rmin - sqrt( dx.^2 + dy.^2 ) );                                % conv. kernel
h = h ./ sum( h( : ) );             % NORMALIZE THE CONVOLUTION KERNEL
dHs = ones( nEl, 1 );
prjMAX = @(v,eta,beta) coth( beta ) .* tanh( v * beta );                   % relaxed MAX step function (limit of the ETA --> 0)
prjMIN = @(v,eta,beta) 1 + coth( beta ) .* tanh( ( v - 1 ) * beta );       % relaxed MIN step function (limit of the ETA --> 1)
prjETA = @(v,eta,beta) ( tanh( beta * eta ) + tanh( beta * ( v( : ) - eta ) ) ) ./ ...
    ( tanh( beta * eta ) + tanh( beta * ( 1 - eta ) ) );                   % relaxed Heaviside projection
dprjMAX = @(v,eta,beta) beta * coth( beta ) .* ( 1 - tanh( v * beta ).^2 );% derivative w.r.t. x (MAX)
dprjMIN = @(v,eta,beta) beta * coth( beta ) .* ( 1 - tanh( ( v - 1 ) * beta ).^2 ); % derivative w.r.t. x (MIN)
dprjETA = @(v,eta,beta) beta * ( 1 - tanh( beta * ( v - eta ) ).^2 ) ./ ...
    ( tanh( beta * eta ) + tanh( beta * ( 1 - eta ) ) );                   % derivative w.r.t. x (ETA)
invETA = @(v,eta,beta) atanh( v * ( tanh( beta * eta ) + ...
    tanh( beta * ( 1 - eta ) ) ) - tanh( beta * eta ) ) / beta + eta;      % inverse ETA projection
switch prjType
    case 'max', prj = prjMAX; dprj = dprjMAX;
    case 'min', prj = prjMIN; dprj = dprjMIN;
    case 'eta', prj = prjETA; dprj = dprjETA;
end
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% LOADS, SUPPORTS AND PASSIVE DOMAINS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% MBB BEAM WITH distributed FORCE (based on 240x80 discretization) >>>>>>>>
szLd                 = nelx / 30;  % number of elements for spreading the load
lcDof                = 2 * nodeNrs( 1, 1 : szLd );
magF                 = 1;
%F                    = fsparse( lcDof( : ), 1, -magF ./ ( szLd - 1 ), [ nDof, 1 ] );
F = sparse( double(lcDof(:)), ones(size(lcDof)), -magF./(szLd-1), nDof,1 );

F( lcDof( 1 ), : )   = 0.5 * F( lcDof( 1 ), : );
F( lcDof( end ), : ) = 0.5 * F( lcDof( end ), : );
% define passive solid and void domains ...................................
[ pasS, pasV ] = deal([],[]);
pasS           = pasS( : );
pasV           = pasV( : );
pasEl          = union( pasS, pasV );
act            = setdiff( ( 1 : nEl )', pasEl );        % set of active DVs
% define displacement boundary conditions .................................
fix  = union( 1 : 2 : 2 * ( nely + 1 ), 2 * nodeNrs( end, end - szLd / 2 + 1 : end ) );
free = setdiff( 1 : nDof, fix );
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% DEFINE ANONYMOUS FUNCTIONS USED IN THE SCRIPT >>>>>>>>>>>>>>>>>>>>>>>>>>>
% ------- function computing the volume fraction for active/passive domains
vfXp = @( v ) sum( v( : ) ) ./ ( length( act ) + length( pasS ) );
% ------- functions used for filtering and back-filtering the density field
applyFilter = @( v, w ) imfilter( reshape( v, nely, nelx ), w, bcF );  % function applying the filter
backFilter  = @( v, w ) imfilter( reshape( v, nely, nelx ), w, bcF );  % function enforcing back-filtering
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% ALLOCATE AND INITIALIZE OTHER PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
[ x, dsK, dsR, dVdxp,  ] = deal( zeros( nEl, 1 ) );                 % initialize vectors of size nEl x 1
% - initialize design variables, density distribution and volume derivative
tmp = (volfrac*(nEl-length(pasV))-length(pasS))./length(act);
% --------------------------------- apply volume fraction on the active set
x( act ) = invETA( tmp, eta, beta );
x( pasS ) = 1;                                                             % set x = 1 on pasS set
dVdxp( act, 1 ) = 1 ./ ( length( act ) + length( pasS ) );                 % xPhys-derivative of volume (constant)
% ----- allocate and initialize other arrays used for tracking optimization
[ xTilde, xPhys ] = deal( x );      % reative and intermediate fields
loop      = 0;                               % initialize iteration counter
clear dx dy;                             % clear some variables from memory
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%% Stochastic model setup
n_disc = 3000;
y = linspace(wInt(1),wInt(2),n_disc);
maxsmpl = 2500;
% mini-batches
btchsz = 4;
if floor(maxsmpl/btchsz) ~= maxsmpl/btchsz
    maxsmpl = floor(maxsmpl/btchsz)*btchsz;
    fprintf('max number of samples was adjusted to %i\n', maxsmpl)
end
x_birth = reshape(repmat(1:maxsmpl/btchsz,btchsz,1),[1, maxsmpl]);
x_ind = 1:maxsmpl;
leavers = 1:btchsz;
X = wInt(1) + rand(1,maxit*btchsz)*(wInt(2)-wInt(1));
y_weight = volfrac*sqrt(nEl);
y_diff = pdist2(y',X');
y_diff = y_diff/max(max(y_diff,1e-10),[],'all');
[Gra,DesH,ComH] = deal( zeros(nEl,maxsmpl), zeros(nEl,maxsmpl), zeros(1,maxsmpl) );
%% ================================================ START OPTIMIZATION LOOP
while loop < maxit
    loop = loop + 1;                                                       % update iteration counter    
    %% COMPUTE PHYSICAL FIELD BY APPLYING FILERING AND PROJECTION >>>>>>>>>
    % ------------------------------- density or sensitivity filtering step
    if ft > 0                                        % apply density filter
        tmp = applyFilter( x, h );
        xTilde( : ) = tmp( : );    % dont' replace xPhys on passive regions
        if prjSel > 0                          % apply nonlinear projection
            dHs = dprj( xTilde, eta, beta );                               % WARNING! this is for the relationship xPhys = f(xTilde)
            xPhys = prj( xTilde, eta, beta );                              % projected (physical) field
        else
            xPhys = xTilde;
        end
    else                                         % apply sensitivity filter
        xTilde = x;
    end
    %% APPLY EXTERNAL PENALIZATION FUNCTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sK0        = intF( xPhys, penal, E0, Emin );   % interpolate stiffness
    dsK( act ) = dintF( xPhys( act ), penal, E0, Emin );
    sR0        = intF( xPhys, penalR, R0, Rmin );     % interpolate density
    dsR( act ) = dintF( xPhys( act ), penalR, R0, Rmin );
    %% PERFORM FINITE ELEMENT ANALYSIS, COMPUTE COMPLIANCE & SENSITIVIY >>>
    % ------------------------------------ assemble global stiffness matrix
    sK = Ke .* sK0( : )';
    %K = fsparse( Iar( :, 1 ), Iar( :, 2 ), sK( : ), [ nDof, nDof ] );      % assemble stiffness matrix (lower half)
    K = sparse( double(Iar(:,1)), double(Iar(:,2)),sK(:), nDof, nDof );
    sM = Me .* sR0( : )';
    %M = fsparse( Iar( :, 1 ), Iar( :, 2 ), sM( : ), [ nDof, nDof ] );      % assemble mass matrix (lower half)
    M = sparse( double(Iar(:,1)), ...
            double(Iar(:,2)), ...
            sM(:), ...
            nDof, nDof );
    % Sample mini-batch
    ulim = min(maxsmpl,loop*btchsz);
    curr_freqs = X(:,(loop-1)*btchsz+1:loop*btchsz);
    g0f = zeros(1,btchsz);
    g1f = zeros(1,btchsz);
    dg0f = zeros(size(x,1),btchsz);
    dg1f = zeros(size(x,1),btchsz);
    parfor ba = 1:btchsz
        U = zeros( nDof, size( F, 2 ) );
        w0 = curr_freqs(:,ba);
        S = K - w0^2 * M;
        S = S + S.'- diag( diag( S ) );
        % ----------------------------------------- solve equilibrium equations
        dS = decomposition( S( free, free ) );                                 % system matrix
        U( free, : ) = dS \ F( free, : );                                      % f/b substitution
        % ---- compliance and compliance sensitivity (for each RHS if multiple)
        % Complex version
        c = abs( F' * U );                        % compute current compliances
        alpha = U'*F / c;
        tmp = dsK .* sum( ( U( cMat ) * Ke0 ) .* U( cMat ), 2 );
        tmp = tmp - dsR .* sum( ( U( cMat ) * Me0 * w0^2 ) .* U( cMat ), 2 );
        dcdxp = - real(alpha*tmp);
        %% DEFINE OBJECTIVE AND CONSTRAINT(S) OF THE OPTIMIZATION PROBLEM >>>>>
        % --- STEP 1) add volume and compliance to the optimization formulation
        csF = volfrac;
        g0  = c;                                                           % objective: maximum compliance
        dg0 = dcdxp;                                                       % objective sensitivity
        g1  = vfXp( xPhys ) / csF - 1;                                     % constraint: maximum volume
        dg1 = dVdxp( : ) / csF;                                            % sensitivity of the constraint
        % ------- STEP 2) apply filtering (and projection) to the sensitivities
        if ft > 0               % backfiltering when the density filter is used
            dg0( : ) = reshape( backFilter( dg0( : ) .* dHs, h ), [], 1 );     % back-project/back-filter objective sensitivity
            for j = 1 : size( dg1, 2 )
                dg1(:,j) = reshape( backFilter(dg1(:,j).*dHs,h), [], 1 );  % back-project/back-filter sensitivity of the constraints
            end
        else                % backfiltering when the sensitivity filter is used
            tmp = reshape( max( 1e-3, x( : ) ), nely, nelx );
            dg0( : ) = reshape( backFilter(x(:).*dg0(:).*dHs,h)./tmp,[],1);
            for j = 1 : size( dg1, 1 )
                dg1(:,j)=reshape(backFilter((dHs.*x(:)).*dg1(:,j),h)./tmp,[],1); % back-project/back-filter sensitivity of the constraints
            end
        end
        %% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        g0f(:,ba) = g0;
        dg0f(:,ba) = dg0;
        g1f(:,ba) = g1;
        dg1f(:,ba) = dg1;
    end
    % store samples
    Gra(:,leavers) = dg0f;
    ComH(leavers) = g0f;
    DesH(:,leavers) = repmat(xPhys,1,btchsz);
    dg1 = dg1f(:,1);
    % CS weights
    [~,csw] = min( vecnorm(xPhys-DesH(:,1:ulim),2,1) + y_weight*y_diff(:,x_ind(1:ulim)), [], 2);
    weights = sum(csw==1:ulim)/n_disc;
    % Determine new leavers
    if (loop+1)*btchsz <= maxsmpl
        leavers = loop*btchsz+1 : (loop+1)*btchsz;
    else
        minweight = sort(weights);
        ind_can = find(weights-minweight(btchsz)<1e-8);
        can_birth = x_birth(ind_can);
        [~,iind] = sort(can_birth);
        leavers = ind_can(iind(1:btchsz));
        x_ind(leavers) = loop*btchsz+1 : (loop+1)*btchsz;
        x_birth(leavers) = loop+1;
    end
    % integrate model
    Compl = sum(weights.*ComH(:,1:ulim),2);Cp = (sum(weights.*(ComH(:,1:ulim).^pnorm),2))^(1/pnorm);
    if mod(loop,25) == 0 || loop == 1
        com0 = Compl;           % adjust normalization
    end
    dc = 1/com0*sum(weights.*(((ComH(:,1:ulim)/com0).^(pnorm-1)).*Gra(:,1:ulim)),2);
    %% UPDATE DESIGN VARIABLES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    xT = x( act );
    [ xU, xL ] = deal( xT + move, xT - move );
    ocP = xT .* real( sqrt( max(1e-10,-dc( act )) ./ dg1(act) ) );
    l = [ 0, (mean( ocP ) / volfrac) ];                                        % initial estimate for LM
    while ( l( 2 ) - l( 1 ) ) / ( l( 2 ) + l( 1 ) ) > 1e-8  && l(2) > 1e-40                 % OC resizing rule
        lmid = 0.5 * ( l( 1 ) + l( 2 ) );
        x( act ) = max( max( min( min( ocP / lmid, xU(act) ), 1 ), xL(act) ), 0 );
        tmp = applyFilter( x, h );
        xTilde( : ) = tmp( : );
        xPhys( act ) = xTilde( act );
        xPhys = prj( xPhys, eta, beta );
        if mean( xPhys ) > volfrac, l( 1 ) = lmid; else, l( 2 ) = lmid; end
    end
    [penal,beta] = deal(cnt(penal,penalCnt,loop), cnt(beta,betaCnt,loop));  % apply conitnuation on parameters
    % ------------------------ RL. 7) PRINT CURRENT RESULTS AND PLOT DESIGN
    fprintf( 'It.:%5i C:%7.4f Cp:%7.4f V:%7.3f penal:%7.2f beta:%7.1f eta:%7.2f \n', ...
        loop, Compl, Cp, mean(xPhys), penal, beta, eta);
    colormap( gray ); imagesc( 1 - reshape( xPhys, nely, nelx ) );
    caxis([0 1]); axis equal off; drawnow
end
end

%% Auxiliary
function [ nodeNrs, nDof, cMat, Ke, Ke0, Me, Me0, Iar, sI, sII ] = discQ4( nely, nelx, nu )

nEl     = nelx * nely;
nodeNrs = int32(reshape(1:(1+nely)*(1+nelx),1+nely,1+nelx));               % node numbering (int32)
lC      = int32( [ 0, 1, 2 * nely + [ 2, 3, 0, 1 ], -2, -1 ] );            % local connectivity
cMat    = reshape( 2 * nodeNrs( 1:end - 1, 1:end - 1 ) + 1, nEl, 1 ) + lC; % connectivity matrix
nDof    = ( 1 + nely ) * ( 1 + nelx ) * 2;                                 % total number of DOFs
% ---------------------------------------------- elemental stiffness matrix
c1 = [12;3;-6;-3;-6;-3;0;3;12;3;0;-3;-6;-3;-6;12;-3;0;-3;-6;3;12;3;...
    -6;3;-6;12;3;-6;-3;12;3;0;12;-3;12];
c2 = [-4;3;-2;9;2;-3;4;-9;-4;-9;4;-3;2;9;-2;-4;-3;4;9;2;3;-4;-9;-2;...
    3;2;-4;3;-2;9;-4;-9;4;-4;-3;-4];
Ke = 1 / ( 1 - nu^2 ) / 24 * ( c1 + nu .* c2 );                            % lower symmetric part of Ke
Ke0( tril( ones( 8 ) ) == 1 ) = Ke';
Ke0 = reshape( Ke0, 8, 8 );
Ke0 = Ke0 + Ke0' - diag( diag( Ke0 ) );                                    % recover full elemental matrix
% --------------------------------------------------- elemental mass matrix
Me = [4,0,2,0,1,0,2,0,4,0,2,0,1,0,2,4,0,2,0,1,0,4,0,2,0,1,4,0,...
    2,0,4,0,2,4,9,4]'./ (36*nelx*nely);                                                % lower symmetric part of Me
Me0( tril( ones( 8 ) ) == 1 ) = Me';
Me0 = reshape( Me0, 8, 8 );
Me0 = Me0 + Me0' - diag( diag( Me0 ) );                                    % recover full mass matrix
% ----------------------- build indices for the assembly of global matrices
[ sI, sII ] = deal( [ ] );
for j = 1 : 8    % build assembly indices for the lower symmetric part of K
    sI = cat( 2, sI, j : 8 );
    sII = cat( 2, sII, repmat( j, 1, 8 - j + 1 ) );
end
[ iK, jK ] = deal( cMat( :, sI )', cMat( :, sII )' );
Iar = sort( [ iK( : ), jK( : ) ], 2, 'descend' );                          % indices for K assembly
%
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