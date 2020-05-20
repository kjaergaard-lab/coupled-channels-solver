function [Eout,results] = MultiChannelBound(InitStateLabel,Ein,Bin,OutputFile,BasisSetFile,opt)
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


%% Fundamental constants and scaling
load(BasisSetFile);

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

scale = const.cm2K*const.K2A(mass); %Converts energies in cm^{-1} to wavenumber^2 in Angstrom^2
Nruns = numel(Bin);
Escale = Ein*const.K2A(mass);  %Convert from K to wavenumber^2 in Angstrom^2
B = Bin*1e-4;                   %Convert from Gauss to Tesla


%% Subspace restriction

if opt.usedipole
    %Restriction to same mL+mF1+mF2=mJ and L'=L-2,L,L+2
    TotalM = sum(InitStateLabel(2:3));
    MatchLabel = [InitStateLabel(1) TotalM];

    [BVintMatch,~] = MatchQuantumNumbers(MatchLabel,[BVint(:,1:2),IntMValues],2:4,'Dipole');
    BVint = BVint(BVintMatch,:);

    [BV1Match,BV1] = MatchQuantumNumbers(MatchLabel,BV1,2:6,'Dipole');
    [BV2Match,BV2] = MatchQuantumNumbers(MatchLabel,BV2,[2,4:6],'Dipole');
    [BV3Match,BV3] = MatchQuantumNumbers(MatchLabel,BV3,[2,4,6],'Dipole');
    [BV4Match,BV4] = MatchQuantumNumbers(MatchLabel,BV4,[2,6],'Dipole');
else
    %Restriction to same mF=mF1+mF2 AND same L AND mL
    MatchLabel = InitStateLabel;

    [BVintMatch,~] = MatchQuantumNumbers(MatchLabel,[BVint(:,1:2),IntMValues],1,2,3:4);
    BVint = BVint(BVintMatch,:);

    [BV1Match,BV1] = MatchQuantumNumbers(MatchLabel,BV1,1,2,3:6);
    [BV2Match,BV2] = MatchQuantumNumbers(MatchLabel,BV2,1,2,4:6);
    [BV3Match,BV3] = MatchQuantumNumbers(MatchLabel,BV3,1,2,[4,6]);
    [BV4Match,BV4] = MatchQuantumNumbers(MatchLabel,BV4,1,2,6); 
end

BT21 = BT21(BV2Match,BV1Match);
BT31 = BT31(BV3Match,BV1Match);
BT43 = BT43(BV4Match,BV3Match);
LMat = diag(BV1(:,1));

NSubChannels = size(BV4,1);
IR = eye(NSubChannels);

%% Hamiltonian settings
%Internal Hamiltonian
tmp1 = 0.5*(BV3(:,3).*(BV3(:,3)+1)-Spin*(Spin+1)-NSpin1*(NSpin1+1));
tmp2 = 0.5*(BV3(:,5).*(BV3(:,5)+1)-Spin*(Spin+1)-NSpin2*(NSpin2+1));
Hint0 = scale*diag(Ahfs1*tmp1+Ahfs2*tmp2);   %In the spin-nucleus coupled basis, in 1/Angstroms^2
Hint0 = BT31'*Hint0*BT31;   %In uncoupled basis, to be added to Zeeman Hamiltonian

%Spin exchange Hamiltonian
SpinProj = diag(0.5*(BV2(:,3).*(BV2(:,3)+1)-2*Spin*(Spin+1)));  %Defined in spin-coupled basis

%Dipole-Dipole Hamiltonian
if opt.usedipole  
    Hdd = -sqrt(6)*FineStructure.^2.*EHartree.*a_Bohr^3/(2*pi*hbar*LightSpeed)*1e-2*DipoleDipole(BV1);    %In uncoupled basis, in cm^{-1}*Angstroms^3
else
    Hdd = zeros(size(SpinProj));
end

%% Create internal transformation matrices
disp('Constructing basis sets');
Hint = zeros(NSubChannels,NSubChannels,Nruns);
BT1int = Hint;
BT2int = Hint;
BScale=mu_B./(2*pi*hbar*LightSpeed)*1e-2;
for nn=1:Nruns
    [Hint(:,:,nn),BT1int(:,:,nn),BT2int(:,:,nn)] = Generate_Hint(Hint0,BV1,BVint,BT21,B(nn),gS_1,gS_2,gI_1,gI_2,Ahfs1,Ahfs2,NSpin1,NSpin2,BScale,scale);
end

%% Prepare all matrix operators in the total S (BV2) basis
Hdd = BT21*Hdd*BT21';
Hint0 = BT21*Hint0*BT21';


%% Loop over magnetic field inputs
results(Nruns,1) = boundresults;
maxBoundStates = 0;
for nn=1:Nruns
    fprintf(1,'Run %d/%d\n',nn,Nruns);
    ops = boundoperators(scale,LMat,SpinProj,Hdd,Hint0,Hint(:,:,nn),BT2int(:,:,nn));
    [r,opt.blocks] = makegrid(@(x) PotentialFunc(x,scale,LMat,SpinProj,Hdd,Hint0)+Hint(:,:,nn),max(Escale),opt);
    Eranges = findranges(r,PotentialFunc,Escale,true,ops,opt);
    if size(Eranges,2)>maxBoundStates
        maxBoundStates = size(Eranges,2);
    end
    
    if opt.debug
        fprintf(1,'Number of bound states between [%#.5g,%#.5g]: %d\n',Ein(1),Ein(2),size(Eranges,2));
    end
    results(nn).r = r;
    results(nn).BT2int = BT2int(:,:,nn);
    results(nn).BV2 = BV2;
    results(nn).BVint = BVint;
    results(nn).basis = 1;
    
    for mm=1:size(Eranges,2)
        if opt.output
            [Ebound,wf] = solvebound(r,PotentialFunc,Eranges(:,mm),ops,opt);
            results(nn).E{mm} = Ebound/const.K2A(mass);
            results(nn).wf{mm,1} = wf;
        else
            Ebound = solvebound(r,PotentialFunc,Eranges(:,mm),ops,opt);
            results(nn).E{mm} = Ebound/const.K2A(mass);
        end
    end
end

Eout = NaN(Nruns,maxBoundStates);
for nn=1:Nruns
    Eout(nn,1:results(nn).count) = cell2mat(results(nn).E);
end

% nn = 1;
% ops = boundoperators(scale,LMat,SpinProj,Hdd,Hint0,Hint(:,:,nn),BT2int(:,:,nn));
% [r,opt.blocks] = makegrid(@(x) PotentialFunc(x,scale,LMat,SpinProj,Hdd,Hint0)+Hint(:,:,nn),max(Escale),opt);

if ~isempty(OutputFile)
    save(OutputFile);
end


end


function [Hint,BT1int,BT2int]=Generate_Hint(Hint0,BV1,BVint,BT21,B,gS_1,gS_2,gI_1,gI_2,Ahfs1,Ahfs2,NSpin1,NSpin2,mu_B,Scale)

NSubChannels=size(Hint0,1);
IR=eye(size(Hint0));

Hint=Hint0+B*mu_B*Scale*diag(gS_1.*BV1(:,3)+gS_2.*BV1(:,5)+gI_1.*BV1(:,4)+gI_2.*BV1(:,6));   %In uncoupled basis

% Determine transformation from basis of internal states in presence of
% B to fully uncoupled basis
[~,V1,~,MagLabel1]=HyperfineSolveInt(B,NSpin1,Ahfs1,gS_1,gI_1,mu_B);  %V is the one particle version of BT1int
V1(log(abs(V1))<-30)=0;
[~,V2,~,MagLabel2]=HyperfineSolveInt(B,NSpin2,Ahfs2,gS_2,gI_2,mu_B);  %V is the one particle version of BT1int
BT1int=zeros(NSubChannels);
for row=1:NSubChannels, 
    for col=1:NSubChannels
        if BV1(row,1)==BVint(col,1)
            tmp1=repmat(BV1(row,:),size(V1,1),1);
            idx1=all(tmp1(:,3:4)==MagLabel1,2);
            tmp2=repmat(BV1(row,:),size(V2,1),1);
            idx2=all(tmp2(:,5:6)==MagLabel2,2);
            BT1int(row,col)=V1(idx1,BVint(col,3))*V2(idx2,BVint(col,4));
        end
    end
end
BT2int=BT21*BT1int;

Hint=BT1int'*Hint*BT1int;   %Internal Hamiltonian in basis of internal states (diagonal)
[~,minIndex] = min(diag(Hint));
Hint=Hint-Hint(minIndex,minIndex)*IR;
        
end

function [E,V,Hint,StateLabel]=HyperfineSolveInt(B,I,Ahfs,gS,gI,mu_B)
S=0.5;   %Electron spin

NI=2*I+1;
NS=2*S+1;
Ndim=NI*NS;
mS=-S:S;
mI=-I:I;
StateLabel=[reshape(repmat(mS,NI,1),Ndim,1) repmat(mI(:),NS,1)];

F=abs(I-S):abs(I+S);
StateLabelF=[];
for nn=1:numel(F)
    mF=(-F(nn):F(nn))';
    StateLabelF=[StateLabelF;[repmat(F(nn),numel(mF),1) mF]];
end

U=zeros(Ndim);
for a=1:Ndim
    F=StateLabelF(a,1);
    mF=StateLabelF(a,2);
    for b=1:Ndim
        mI=StateLabel(b,2);
        mS=StateLabel(b,1);
        if abs(mI+mS)>F || (mI+mS)~=mF
            continue;
        else
            U(a,b)=ClebschGordan(I,S,F,mI,mS,mF);
        end
    end
end


%% Bare Hamiltonian
H=zeros(Ndim);
F=abs(I-S):abs(I+S);
for a=1:Ndim
    mI1=StateLabel(a,2);
    mS1=StateLabel(a,1);
    for b=1:Ndim
        mI2=StateLabel(b,2);
        mS2=StateLabel(b,1);
        for c=1:numel(F)
            if abs(mI1+mS1)>F(c) || abs(mI2+mS2)>F(c) || (mI1+mS1)~=(mI2+mS2)
                continue;
            end
                CG=ClebschGordan(I,S,F(c),mI1,mS1,mI1+mS1).*ClebschGordan(I,S,F(c),mI2,mS2,mI2+mS2);
                K=F(c)*(F(c)+1)-I*(I+1)-S*(S+1);
                tmp=CG*Ahfs/2*K;
                H(a,b)=H(a,b)+tmp;
        end
    end
end



%% Zeeman field
HB=diag(mu_B.*B*(gS*StateLabel(:,1)+gI*StateLabel(:,2)));
Hint=H+HB;
[V,E]=eig(Hint);
D=diag(E);
[D,idx]=sort(D);
V=V(:,idx);
E=diag(D);



end
