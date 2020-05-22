function [Output,wf,S] = MultiChannel(initLabel,Ein,Bin,outputFile,basis,opt)
% MultiChannel computes the scattering properties of a pair of alkali metal
% atoms
%   Usage 1: MultiChannel(InitStateLabel,Ein,Bin,OutputFile,BasisSetFile,DipoleFlag,IntParams)
%
%   InitStateLabel is four-element row vector [L,mL,State1,State2]
%   specifying the entrance channel.  L is the angular momentum, mL is its
%   projection, State1 and State2 are the internal states of the atoms,
%   labelled in increasing ground state energy in the presence
%   of a magnetic field.  For identical particles, specify the lowest
%   number first.  For KRb collisions, State1 is for K and State2 is for Rb
%
%   Ein and Bin specify entrance channel energies and magnetic fields in uK
%   and G, respectively.  Only one of these can have more than one element
%
%   OutputFile is the name of the output .mat file where the results are
%   saved.  If empty, the results are not saved
%
%   BasisSetFile is the name of the .mat file that specifies the species,
%   constants, potential functions, and basis vectors and their
%   transformations
%
%   DipoleFlag=0,1 sets whether or not to include the magnetic
%   dipole-dipole interaction between valence electrons.  Including it is
%   more accurate but uses more resources.
%
%   IntParams is a structure that sets integration parameters.  rMin is the
%   starting distance of integration in Angstroms (cannot be 0). rMax is
%   the final distance - 500 angstroms seems to work well. dr_scale sets
%   the spatial step size as a fraction of the WKB wavelength.  Smaller
%   numbers are more accurate but take longer. drMin and drMax set minimum
%   and maximum step sizes. ParSet=0,1 sets whether or not to use parallel
%   processing.  Default values are
%   IntParams.rMin=1e-1;
%   IntParams.rMax=500;
%   IntParams.dr_scale=1e-2;
%   IntParams.drMax=500;
%   IntParams.drMin=1e-3;
%   IntParams.ParSet=1;
%
%   Usage 2: Output=MultiChannel(InitStateLabel,Ein,Bin,OutputFile,BasisSetFile,DipoleFlag,IntParams)
%   Output is a structure with fields: k, S, T, S2, T2, and OutIdx
%
%   See also BasisSetGeneration, CalculatePartialCrossSection, FitResonance


%% Load basis set and restrict to particular subspace
if isa(basis,'atompair')
    basis = copy(basis);
elseif ischar(basis)
    V = load(basisSetFile,'basis');
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
T = S-repmat(eye(basis.Nchannels),[1,1,Nruns]); %T-matrix is just S-1
if basis.symmetry~=0
    [Ssym,BVSymInt,ia] = SymmetrizeSMatrixArbL(S,basis.bvint,basis.symmetry);
    if opt.getwf
        for nn=1:numel(wf)
            wf{nn}.u = wf{nn}.u(ia,:);
        end
    end
    tmp = size(Ssym);
    if Nruns==1
        tmp(3) = 1;
    end
    Tsym = Ssym-repmat(eye(tmp(1:2)),[1,1,tmp(3)]);
    OutIdx = find(all(BVSymInt(:,[1,2,5,6])==repmat([initLabel(1:2),sort(initLabel(3:4),2)],tmp(1),1),2));
    
    if Nruns==1
        S2 = Ssym(OutIdx,:).';
        T2 = Tsym(OutIdx,:).';
    else
        S2 = shiftdim(Ssym,2);
        S2 = squeeze(S2(:,:,OutIdx));
        T2 = shiftdim(Tsym,2);
        T2 = squeeze(T2(:,:,OutIdx));
    end
    Output.bvsym = BVSymInt;
    Output.S = Ssym;
    Output.T = Tsym;
else    
    OutIdx = basis.findstate(basis.bvint(:,[1,2,5,6]),initLabel);
    
    if Nruns==1
        S2 = S(OutIdx,:).';
        T2 = T(OutIdx,:).';
    else
        S2 = shiftdim(S,2);
        S2 = squeeze(S2(:,:,OutIdx));
        T2 = shiftdim(T,2);
        T2 = squeeze(T2(:,:,OutIdx));
    end
    Output.bvint = basis.bvint;
    Output.S = S;
    Output.T = T;
end

if eflag
    k = sqrt(2*basis.mass/const.hbar^2*E(:)*1e-6*const.kb);
else
    k = sqrt(2*basis.mass/const.hbar^2*E(:)*1e-6*const.kb);
end

Output.k = k;
Output.S2 = S2;
Output.T2 = T2;
Output.OutIdx = OutIdx;
Output.mass = basis.mass;
Output.Ein = Ein;
Output.Bin = Bin;

if size(S,3)>1
    sm = ScatteringMatrix(S,basis.bvint,basis.symmetry,initLabel);
    sm.mass = basis.mass;
    sm.E = Ein;
    sm.B = Bin;

    Output = sm;
end

if ~isempty(outputFile)
    save(outputFile);
end

end

