function mout = manolopoulos_scatt(E,ops,basis,opt)
%manolopoulos_scatt Integrates the matrix Riccati equation using the
%improved log-derivative method of D. E. Manolopoulos for scattering
%problems
%
%   mout = manolopoulos(E,ops,basis,opt) integrates the matrix Riccati
%   equation outwards assuming an energy E.
%
%   ops is an instance of ALKALIOPERATORS
%
%   basis is an instance of ATOMPAIR
%
%   opt is an instance of SCATTOPTIONS
%
%   mout is either the scattering matrix, when opt.getwf = false, or a
%   structure with fields S (scattering matrix) and wf (radial wavefunction
%   structures).  wf is a structure with fields r (radial position) and u
%   (radial wavefunctions for each channel)

%% Prep
I = eye(ops.Nch);   %Define the identity
%
% Set the initial value of the log-derivative matrix to be diagonal with
% values equal to the WKB wave numbers
%
Ynew = sqrt(ops.potential(opt.rmin,E)).*I;
if opt.getwf
    %
    % If the user wants the wavefunctions, certain propagation matrices
    % need to be stored
    %
    Zout = zeros(ops.Nch);
    rout = opt.rmin;
end

%% Calculate adiabatic potentials
%
% The adiabatic potentials are used to determine the step size
%
ra = linspace(opt.rmin,opt.rmax,5e3);
f = ops.potential(ra,E);
adiabat = zeros(ops.Nch,numel(ra));
for kk = 1:numel(ra)
    adiabat(:,kk) = sort(eig(f(:,:,kk)));
end

%% Perform integration
%
% The "change flag" indicates when the basis switches from the spin-coupled
% basis suitable for the BO potentials to the "internal" basis suitable for
% the asymptotic states.
%
% Integration is done in blocks to take advantage of pre-allocation
%
changeFlag = false;
numBlocks = floor(opt.rmax/opt.blocksize);

for jj = 1:numBlocks
    %
    % Get start and end radial positions
    %
    rstart = opt.rmin+(jj-1)*opt.blocksize;
    rend = rstart+opt.blocksize;
    %
    % Calculate the maximum wavenumber from the lowest adiabat, and update
    % the step size accordingly. If all channels are closed in this block,
    % then this will keep the last value of kmax
    %
    kmax = adiabat(1,ra>=rstart & ra<=rend);
    kmax = kmax(kmax<0);
    kmax = max(sqrt(abs(kmax)));
    if ~isempty(kmax)
        dr = opt.drscale*2*pi./kmax;
        dr = max(min(dr,opt.drmax),opt.drmin);
        dr = min(dr,opt.blocksize/2);
    end
    %
    % Get number of points and define a position grid with half-step size
    % h. Pre-allocate the necessary arrays
    %
    N = max(ceil((rend-rstart)/dr),2);
    rb = linspace(rstart,rend,N);
    h = (rb(2)-rb(1))/2;
    
    if opt.getwf
        Z = zeros(ops.Nch,ops.Nch,N);   %Propagator matrices
    end
    %
    % Potential matrices at position and half-step advanced position
    %
    M = ops.potential(rb,E);
    M2 = ops.potential(rb+h,E);
    
    %% Invariant embedding solution
    for nn = 1:numel(rb)-1
        %
        % See D. E. Manolopolous's paper J. Chem. Phys. 85, 6425 (1986); https://doi.org/10.1063/1.451472
        % for the original source. See also A. E. Thornley's paper J. Chem. Phys. 101, 5578 (1994); https://doi.org/10.1063/1.467345
        % for how to get the wave functions
        %
        p2 = diag(M2(:,:,nn));
        p = sqrt(abs(p2));
        y1 = p.*coth(p*h).*(p2>0)+p.*cot(p*h).*(p2<0)+1./h.*(p2==0);
        y2 = p.*csch(p*h).*(p2>0)+p.*csc(p*h).*(p2<0)+1./h.*(p2==0);
        
        y1 = diag(y1);
        y2 = diag(y2);
        Mref = diag(p2);
        
        Qa = h/3*(M(:,:,nn)-Mref);
        Qc = 4/h*((I-h^2/6*(M2(:,:,nn)-Mref))\I)-4/h*I;
        Qb = h/3*(M(:,:,nn+1)-Mref);
        
        Zac = (Ynew+y1+Qa)\y2;
        Yc = (y1+Qc)-y2*Zac;
        Zcb = (Yc+y1+Qc)\y2;
        Ynew = (y1+Qb)-y2*Zcb;
        
        if opt.getwf
            %
            % This Z matrix is used to back-propagate the wavefunctions
            %
            Z(:,:,nn+1) = Zac*Zcb;
        end
    end
    
    if opt.getwf
        %
        % Assign back-propagators and positions in blocks
        %
        Zout(:,:,end+1:end+size(Z,3)-1) = Z(:,:,2:end);
        rout = [rout rb(2:end)]; %#ok<AGROW>
    end
    
    if ~changeFlag && rend>=10
        %
        % If in the spin-coupled basis (changeFlag = false) and we're
        % farther than 10 Angstroms, change the basis to "internal"
        changeFlag = true;
        ops = ops.rotate;
        Ynew = ops.U*Ynew*ops.U';
        if opt.getwf
            %
            % Note that the back-propagators have to be rotated to the new
            % basis
            %
            for kk=1:size(Zout,3)
                Zout(:,:,kk) = ops.U*Zout(:,:,kk)*ops.U';
            end
        end
    end
end

%% Calculate S-matrix
b = rb(end);                            %Last position
Eabs = diag((E+ops.refE)*I-ops.Hint);   %Energy relative to the open channel
idx = Eabs>0;                           %Determine which channels are open
Eabs = diag(Eabs);
k = sqrt(Eabs(idx,idx));                %Wave numbers for OPEN channels
L = ops.L(idx,idx);
Ynew = Ynew(idx,idx);
%
% This is the definition of the S matrix for open channels. See Ryan
% Thomas's PhD thesis, University of Otago, Eq. 2.62
%
Stmp = (k*RiccatiBessel(b,k,L,2,1)-RiccatiBessel(b,k,L,2,0)*Ynew)/(k*RiccatiBessel(b,k,L,1,1)-RiccatiBessel(b,k,L,1,0)*Ynew);
%
% Technically the S-matrix is undefined for closed channels, but to ensure
% that as the energy increases and new channels become open that the S
% matrix has the same dimension, define the closed elastic channels to have
% an S-matrix element of 1.  This means that the T-matrix for these
% channels is 0, and thus they experience no scattering
%
S = I;
S(idx,idx) = Stmp;

%% Calculate wavefunctions
if opt.getwf
    L = ops.L;
    k = sqrt(Eabs);
    u = zeros(ops.Nch,size(Zout,3));
    %
    % Calculate the "symmetrising" matrix U, which should account for
    % bosonic/fermionic symmetry in the wave-functions.
    %
    if basis.symmetry == 0
        U = I;
    else
        U = I;
        for row = 1:ops.Nch
            for col = 1:ops.Nch
                if all(basis.bvint(row,1:2) == basis.bvint(col,1:2),2)
                    if all(basis.bvint(row,5:6) == basis.bvint(col,5:6),2)
                        U(row,col) = U(row,col) + (2*(1+(basis.bvint(row,5)==basis.bvint(row,6))))^-0.5;
                    end

                    if all(basis.bvint(row,5:6) == basis.bvint(col,[6,5]),2)
                        U(row,col) = U(row,col) + (-1)^basis.bvint(row,1)*(2*(1+(basis.bvint(row,5)==basis.bvint(row,6))))^-0.5;
                    end
                end
            end
        end
    end
    %
    % Phi is the matrix of solutions in the basis of spherical Hankel
    % functions.  It should be zero for closed channels and related as
    % below for open channels defined by IDX.  See R. Thomas, PhD thesis,
    % Eq. 2.59
    %
    phi = zeros(ops.Nch);
    phi(idx,idx) = RiccatiBessel(b,k(idx,idx),L(idx,idx),2,0)-RiccatiBessel(b,k(idx,idx),L(idx,idx),1,0)*Stmp;
    %
    % This is the vector of coefficients that projects out the correct
    % superposition of possible solutions PHI.  See R. Thomas, PhD thesis,
    % Eq. 2.66
    %
    C = U*diag(-sqrt(pi*(2*L+I)/k)*(1i)^(L-I));
    %
    % Back-propagate the wavefunctions
    %
    u(:,end) = phi*C;
    for nn=size(Zout,3):-1:2
        phi = Zout(:,:,nn)*phi;
        u(:,nn-1) = phi*C;
    end
    wf.r = rout;
    wf.u = u;
end

if opt.getwf
    mout.S = S;
    mout.wf = wf;
else
    mout = S;
end

end
    
