function [sm,wf] = MultiChannel(initLabel,Ein,Bin,outputFile,basis,opt)
% MultiChannel Computes the scattering properties of a pair of alkali metal
% atoms
%   Usage 1: MultiChannel(initLabel,Ein,Bin,outputFile,basis,opt)
%
%   initLabel is four-element row vector [L,mL,State1,State2]
%   specifying the entrance channel.  L is the angular momentum, mL is its
%   projection, State1 and State2 are the internal states of the atoms,
%   labelled in increasing ground state energy in the presence
%   of a magnetic field.  For identical particles, specify the lowest
%   number first.  For KRb collisions, State1 is for K and State2 is for Rb
%
%   Ein and Bin specify entrance channel energies and magnetic fields in uK
%   and G, respectively.  Only one of these can have more than one element
%
%   outputFile is the name of the output .mat file where the results are
%   saved.  If empty, the results are not saved
%
%   basis is either the file name of the .mat file storing the basis
%   information or is a variable of type ATOMPAIR containing the basis
%   information
%
%   opt is an instance of SCATTOPTIONS containing options for the
%   integration of the coupled channels equations
%
%   Usage 2: [sm,wf] = MultiChannel(initLabel,Ein,Bin,outputFile,basis,opt)
%   sm is an instance of the ScatteringMatrix class which contains all the
%   scattering information.  wf are the output radial wavefunctions which
%   are calculated when opt.getwf = true.  It an array of structures with
%   fields r (radial position) and u (wavefunctions) with dimensions 
%   (Nchannels x numel(r)) 
%


%% Load basis set and restrict to particular subspace
if isa(basis,'atompair')
    basis = copy(basis);
elseif ischar(basis)
    V = load(basis,'basis');
    basis = V.basis;
else
    error('Basis set is neither a file nor a variable of type ''atompair''');
end
basis.restrict(initLabel,opt.dipole);

%% Scale input energies and magnetic fields
if any(Bin<1e-3)
    warning('Minimum Bin is 1 mG.  Coercing values...');
    Bin(Bin<1e-3)=1e-3;
end

if any(Ein<=0)
    error('Energies must be positive!');
end

if numel(Ein)>1 && numel(Bin)==1
    Nruns = numel(Ein);
    eflag = true;
elseif numel(Ein)==1 && numel(Bin)>=1
    Nruns = numel(Bin);
    eflag = false;
else
    error('Only one of Bin and Ein can have more than one element!');
end
E = Ein*1e-6*const.K2A(basis.mass); %Convert from uK to inverse wavenumbers
B = Bin*1e-4;                       %Convert from G to T

%% Generate operators/terms in Hamiltonian
disp('Constructing basis sets');
if eflag
    ops = basis.makeOperators(B,initLabel,opt.dipole);
else
    for nn=1:numel(B)
        ops(nn,1) = basis.makeOperators(B(nn),initLabel,opt.dipole); %#ok<AGROW>
    end
end

%% Loop over inputs
disp('Basis sets constructed.  Starting integration');
S = zeros(basis.Nchannels,basis.Nchannels,Nruns);
if opt.parallel
    parfor mm=1:Nruns
        if eflag
            mout = manolopoulos(E(mm),ops,basis,opt);
        else
            mout = manolopoulos(E,ops(mm),basis,opt); %#ok<PFBNS>
        end
        
        if opt.getwf            
            S(:,:,mm) = mout.S;
            wf(mm,1) = mout.wf;
        else
            S(:,:,mm) = mout;
        end
    end
else
    for mm=1:Nruns
        if eflag
            mout = manolopoulos(E(mm),ops,basis,opt);
        else
            mout = manolopoulos(E,ops(mm),basis,opt);
        end
        
        if opt.getwf            
            S(:,:,mm) = mout.S;
            wf(mm,1) = mout.wf; %#ok<AGROW>
        else
            S(:,:,mm) = mout;
        end
    end
end

%% Calculate output quantities
if basis.symmetry~=0
    [Ssym,BVSymInt,ia] = SymmetrizeSMatrixArbL(S,basis.bvint,basis.symmetry);
    if opt.getwf
        for nn=1:numel(wf)
            wf(nn).u = wf(nn).u(ia,:);
        end
    end
    tmp = size(Ssym);
    if Nruns==1
        tmp(3) = 1;
    end
    OutIdx = find(all(BVSymInt(:,[1,2,5,6])==repmat([initLabel(1:2),sort(initLabel(3:4),2)],tmp(1),1),2));
else    
    OutIdx = basis.findstate(basis.bvint(:,[1,2,5,6]),initLabel);

end

if basis.symmetry==0
    sm = ScatteringMatrix(S,basis.bvint,basis.symmetry,OutIdx);
else
    sm = ScatteringMatrix(Ssym,BVSymInt,basis.symmetry,OutIdx);
end
sm.mass = basis.mass;
sm.E = Ein;
sm.B = Bin;

if ~opt.getwf
    wf = [];
end

if ~isempty(outputFile)
    save(outputFile);
end

end

