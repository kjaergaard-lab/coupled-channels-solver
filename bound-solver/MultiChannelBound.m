function [Eout,results] = MultiChannelBound(initLabel,Ein,Bin,outputFile,basis,opt)
% MultiChannelBound Computes the bound state properties of a pair of alkali
% metal atoms
%   Usage 1: MultiChannelBound(initLabel,Ein,Bin,outputFile,basis,opt)
%
%   initLabel is four-element row vector [L,mL,mF]
%   specifying the entrance channel.  L is the angular momentum, mL is its
%   projection, and mF is the total spin projection (mF1+mF2).
%
%   Ein and Bin specify entrance channel energies and magnetic fields in uK
%   and G, respectively.  Only one of these can have more than one element
%
%   outputFile is the name of the output .mat file where the results are
%   saved.  If empty, the results are not saved
%
%   basis is either the file name of the .mat file storing the basis
%   information or is a variable of type ATOMPAIRBASIS containing the basis
%   information
%
%   opt is an instance of BOUNDOPTIONS containing options for the
%   integration of the coupled channels equations
%
%   Usage 2: [sm,wf] = MultiChannelBound(initLabel,Ein,Bin,outputFile,basis,opt)
%   sm is an instance of the ScatteringMatrix class which contains all the
%   scattering information.  wf are the output radial wavefunctions which
%   are calculated when opt.getwf = true.  It an array of structures with
%   fields r (radial position) and u (wavefunctions) with dimensions 
%   (Nchannels x numel(r)) 
%


%% Load basis set and restrict to particular subspace
if isa(basis,'atompairbasis')
    basis = copy(basis);
elseif ischar(basis)
    V = load(basis,'basis');
    basis = V.basis;
else
    error('Basis set is neither a file nor a variable of type ''atompairbasis''');
end
basis.restrict(initLabel,opt.dipole);

%% Scale magnetic fields and energies
if any(Bin<1e-3)
    warning('Minimum Bin is 1 mG.  Coercing values...');
    Bin(Bin<1e-3) = 1e-3;
end

if any(Ein>0)
    error('Energies must be non-positive!');
elseif numel(Ein) ~= 2
    error('Ein must be a two-element vector bracketing the energy range in which to search for bound states');
end

if ~isa(opt,'boundoptions')
    error('Opt must of type boundoptions');
end

Nruns = numel(Bin);
E = Ein*1e-6*scattconst.K2A(basis.mass); %Convert from uK to inverse wavenumbers
B = Bin*1e-4;                       %Convert from G to T

%% Generate operators/terms in Hamiltonian
disp('Constructing basis sets');
for nn=1:numel(B)
    ops(nn,1) = basis.makeOperators(B(nn),initLabel,opt.dipole,'bound'); %#ok<AGROW>
end

opt = repmat(opt,Nruns,1);

%% Loop over magnetic field inputs
disp('Basis sets constructed.  Starting integration');
results(Nruns,1) = boundresults;
% maxBoundStates = 0;
if opt(1).parallel
    parfor nn=1:Nruns
        fprintf(1,'Run %d/%d at B = %.5e T\n',nn,Nruns,B(nn));
        results(nn) = getResults(E,basis,ops(nn),opt(nn));
    end
else
    for nn=1:Nruns
        fprintf(1,'Run %d/%d at B = %.5e T\n',nn,Nruns,B(nn));
        results(nn) = getResults(E,basis,ops(nn),opt(nn));
    end
end

maxBoundStates = 0;
for nn=1:Nruns
    if results(nn).count()>maxBoundStates
        maxBoundStates = results(nn).count();
    end
end

Eout = NaN(Nruns,maxBoundStates);
for nn=1:Nruns
    Eout(nn,1:results(nn).count) = cell2mat(results(nn).E);
end

if ~isempty(outputFile)
    save(outputFile);
end


end

function results = getResults(E,basis,ops,opt)

results = boundresults;
[r,opt.blocks] = makegrid(max(E),ops,opt);
Eranges = findranges(r,E,true,ops,opt);

if opt.debug
    Ein = E/(1e-6*scattconst.K2A(basis.mass));
    fprintf(1,'Number of bound states between [%#.5g,%#.5g]: %d\n',Ein(1),Ein(2),size(Eranges,2));
end
results.r = r;
results.bt2int = ops.BT2int;
results.bv2 = basis.bv2;
results.bvint = basis.bvint;
results.basis = 1;

for mm=1:size(Eranges,2)
    if opt.output
        [Ebound,wf,~,nodes] = solvebound(r,Eranges(:,mm),ops,opt);
        results.E{mm} = Ebound/scattconst.K2A(basis.mass)*1e6;
        results.wf{mm,1} = wf;
        results.nodes{mm} = nodes;
    else
        Ebound = solvebound(r,Eranges(:,mm),ops,opt);
        results.E{mm} = Ebound/scattconst.K2A(basis.mass)*1e6;
    end
end


end

