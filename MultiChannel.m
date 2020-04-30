function [Output,wf]=MultiChannel(InitStateLabel,Ein,Bin,OutputFile,BasisSetFile,DipoleFlag,IntParams)
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
    Bin(Bin<1e-3)=1e-3;
end

if any(Ein<=0)
    error('Energies must be positive!');
end

if numel(Ein)>1 && numel(Bin)==1
    Nruns=numel(Ein);
    EnergyVec=Ein*kb*1e-6;
    EVecScale=1e-20*2*mass/hbar^2*EnergyVec;
    B=Bin*1e-4;
    EFlag=1;
elseif numel(Ein)==1 && numel(Bin)>=1
    Nruns=numel(Bin);
    BVec=Bin*1e-4;
    Energy=Ein*kb*1e-6;
    EScale=1e-20*2*mass/hbar^2*Energy;
    EFlag=0;
else
    error('Only one of Bin and Ein can have more than one element!');
end
Scale=1e-20*1e2*4*pi*mass/hbar*LightSpeed;  %Converts energies given in cm^{-1} to square wavenumbers given in 1/Angstroms^2

%% Integration parameters default settings
if ~isfield(IntParams,'rMin')
    IntParams.rMin=1e-1;
end
if ~isfield(IntParams,'rMax')
    IntParams.rMax=500;
end
if ~isfield(IntParams,'dr_scale')
    IntParams.dr_scale=1e-2;
end
if ~isfield(IntParams,'drMax')
    IntParams.drMax=500;
end
if ~isfield(IntParams,'drMin')
    IntParams.drMin=1e-3;
end
if ~isfield(IntParams,'ParSet')
    IntParams.ParSet=1;
end
if ~isfield(IntParams,'getWavefunctions')
    IntParams.getWavefunctions = false;
end

%% Subspace restriction
InitStateInt=FindState(BVint,InitStateLabel);

if DipoleFlag
    %Restriction to same mL+mF1+mF2=mJ and L'=L-2,L,L+2
    TotalM=InitStateLabel(2)+sum(IntMValues(InitStateInt,:),2);
    MatchLabel=[InitStateLabel(1) TotalM];

    [BVintMatch,~]=MatchQuantumNumbers(MatchLabel,[BVint(:,1:2),IntMValues],2:4,'Dipole');
    BVint=BVint(BVintMatch,:);
    InitStateInt=FindState(BVint,InitStateLabel);

    [BV1Match,BV1]=MatchQuantumNumbers(MatchLabel,BV1,2:6,'Dipole');
    [BV2Match,BV2]=MatchQuantumNumbers(MatchLabel,BV2,[2,4:6],'Dipole');
    [BV3Match,BV3]=MatchQuantumNumbers(MatchLabel,BV3,[2,4,6],'Dipole');
    [BV4Match,BV4]=MatchQuantumNumbers(MatchLabel,BV4,[2,6],'Dipole');
else
    %Restriction to same mF=mF1+mF2 AND same L AND mL
    TotalM=sum(IntMValues(InitStateInt,:),2);

    MatchLabel=[InitStateLabel(1:2) TotalM];

    [BVintMatch,~]=MatchQuantumNumbers(MatchLabel,[BVint(:,1:2),IntMValues],1,2,3:4);
    BVint=BVint(BVintMatch,:);
    InitStateInt=FindState(BVint,InitStateLabel);

    [BV1Match,BV1]=MatchQuantumNumbers(MatchLabel,BV1,1,2,3:6);
    [BV2Match,BV2]=MatchQuantumNumbers(MatchLabel,BV2,1,2,4:6);
    [BV3Match,BV3]=MatchQuantumNumbers(MatchLabel,BV3,1,2,[4,6]);
    [BV4Match,BV4]=MatchQuantumNumbers(MatchLabel,BV4,1,2,6); 
end

BT21=BT21(BV2Match,BV1Match);
BT31=BT31(BV3Match,BV1Match);
BT43=BT43(BV4Match,BV3Match);
LMat=diag(BV1(:,1));

NSubChannels=size(BV4,1);
IR=eye(NSubChannels);

%% Hamiltonian settings
%Internal Hamiltonian
tmp1=0.5*(BV3(:,3).*(BV3(:,3)+1)-Spin*(Spin+1)-NSpin1*(NSpin1+1));
tmp2=0.5*(BV3(:,5).*(BV3(:,5)+1)-Spin*(Spin+1)-NSpin2*(NSpin2+1));
Hint0=Scale*diag(Ahfs1*tmp1+Ahfs2*tmp2);   %In the spin-nucleus coupled basis, in 1/Angstroms^2
Hint0=BT31'*Hint0*BT31;   %In uncoupled basis, to be added to Zeeman Hamiltonian

%Spin exchange Hamiltonian
SpinEx=diag(0.5*(BV2(:,3).*(BV2(:,3)+1)-2*Spin*(Spin+1)));  %Defined in spin-coupled basis

%Dipole-Dipole Hamiltonian
if DipoleFlag 
%     Hdd=-mu0/(4*pi)*(gS*mu_B)^2/(2*pi*hbar*LightSpeed)*1e-2*1e30*DipoleDipole(BV1);    %In uncoupled basis, in cm^{-1}*Angstroms^3
    Hdd=-sqrt(6)*FineStructure.^2.*EHartree.*a_Bohr^3/(2*pi*hbar*LightSpeed)*1e-2*DipoleDipole(BV1);
else
    Hdd=zeros(size(SpinEx));
end

%% Create internal transformation matrices
disp('Constructing basis sets');
Hint=zeros(NSubChannels,NSubChannels,Nruns);
BT1int=Hint;
BT2int=Hint;
E=Hint;
BScale=mu_B./(2*pi*hbar*LightSpeed)*1e-2;
for nn=1:Nruns
    if ~EFlag
        [Hint(:,:,nn),BT1int(:,:,nn),BT2int(:,:,nn)]=Generate_Hint(Hint0,BV1,BVint,BT21,BVec(nn),gS_1,gS_2,gI_1,gI_2,Ahfs1,Ahfs2,NSpin1,NSpin2,BScale,Scale,InitStateInt);
        E(:,:,nn)=EScale*IR-Hint(:,:,nn);
    else
        [Hint(:,:,nn),BT1int(:,:,nn),BT2int(:,:,nn)]=Generate_Hint(Hint0,BV1,BVint,BT21,B,gS_1,gS_2,gI_1,gI_2,Ahfs1,Ahfs2,NSpin1,NSpin2,BScale,Scale,InitStateInt);
        E(:,:,nn)=EVecScale(nn)*IR-Hint(:,:,nn);
    end
end

S=repmat(IR,[1,1,Nruns]);
if DipoleFlag
    BlockSize=5;       %Size of integration blocks in Angstroms
else
    BlockSize=15;      %Size of integration blocks in Angstroms
end

%% Loop over inputs
disp('Basis sets constructed.  Starting integration');
getWavefunctions = IntParams.getWavefunctions;
if IntParams.ParSet==1
    parfor mm=1:Nruns
        if getWavefunctions
            [S(:,:,mm),wf{mm}]=ManolopoulosIntegration(LMat,SpinEx,Hdd,Hint0,E(:,:,mm),BT1int(:,:,mm),BT2int(:,:,mm),IntParams,BlockSize,Scale,BasisSetFile);
        else
            S(:,:,mm)=ManolopoulosIntegration(LMat,SpinEx,Hdd,Hint0,E(:,:,mm),BT1int(:,:,mm),BT2int(:,:,mm),IntParams,BlockSize,Scale,BasisSetFile);
        end
    end
else
    for mm=1:Nruns
        if getWavefunctions
            [S(:,:,mm),wf{mm}]=ManolopoulosIntegration(LMat,SpinEx,Hdd,Hint0,E(:,:,mm),BT1int(:,:,mm),BT2int(:,:,mm),IntParams,BlockSize,Scale,BasisSetFile);
        else
            S(:,:,mm)=ManolopoulosIntegration(LMat,SpinEx,Hdd,Hint0,E(:,:,mm),BT1int(:,:,mm),BT2int(:,:,mm),IntParams,BlockSize,Scale,BasisSetFile);
        end
    end
end

%% Calculate output quantities
Output.BVint = BVint;
T=S-repmat(IR,[1,1,Nruns]); %T-matrix is just S-1
if IdenticalParticles~=0
    [Ssym,BVSymInt,ia]=SymmetrizeSMatrixArbL(S,BVint,IdenticalParticles);
    if IntParams.getWavefunctions
        for nn=1:numel(wf)
            wf{nn}.u = wf{nn}.u(ia,:);
        end
    end
    tmp=size(Ssym);
    if Nruns==1
        tmp(3) = 1;
    end
    Tsym=Ssym-repmat(eye(tmp(1:2)),[1,1,tmp(3)]);
    OutIdx=find(all(BVSymInt==repmat([InitStateLabel(1:2),sort(InitStateLabel(3:4),2)],tmp(1),1),2));
    if EFlag
        k=sqrt(2*mass/hbar^2*EnergyVec(:));   
    else
        k=sqrt(2*mass/hbar^2*Energy(:));  
    end
    if Nruns==1
        S2 = Ssym(OutIdx,:).';
        T2 = Tsym(OutIdx,:).';
    else
        S2=shiftdim(Ssym,2);
        S2=squeeze(S2(:,:,OutIdx));
        T2=shiftdim(Tsym,2);
        T2=squeeze(T2(:,:,OutIdx));
    end
    Output.BVSymInt=BVSymInt;
    Output.S=Ssym;
    Output.T=Tsym;
else    
    if EFlag
        k=sqrt(2*mass/hbar^2*EnergyVec(:));            
    else
        k=sqrt(2*mass/hbar^2*Energy(:));  
    end
    OutIdx=find(all(BVint==repmat(InitStateLabel,NSubChannels,1),2));
    if Nruns==1
        S2 = S(OutIdx,:).';
        T2 = T(OutIdx,:).';
    else
        S2=shiftdim(S,2);
        S2=squeeze(S2(:,:,OutIdx));
        T2=shiftdim(T,2);
        T2=squeeze(T2(:,:,OutIdx));
    end
%     Output.BVint=BVint;
    Output.S=S;
    Output.T=T;
end

Output.k=k;
Output.S2=S2;
Output.T2=T2;
Output.OutIdx=OutIdx;
Output.mass=mass;
Output.Ein=Ein;
Output.Bin=Bin;

if size(S,3)>1
    sm = ScatteringMatrix(S,BVint,IdenticalParticles,InitStateLabel);
    sm.mass = mass;
    sm.E = Ein;
    sm.B = Bin;

    Output = sm;
end

if numel(OutputFile)>0
    save(OutputFile);
end

end

function [S,wf]=ManolopoulosIntegration(LMat,SpinEx,Hdd,H0,E,BT1int,BT2int,IntParams,BlockSize,Scale,BasisFile)

%% Preamble
rMin=IntParams.rMin;
rMax=IntParams.rMax;
IR=eye(size(SpinEx));
NSubChannels=size(IR,1);
load(BasisFile,'PotentialFunc');

%% Define operators in internal basis
SpinExInt=BT2int'*SpinEx*BT2int;  %Spin exchange in coupled total spin basis
HddInt=BT1int'*Hdd*BT1int;    %In spin-coupled basis
%E is already in internal basis
H0Int=BT1int'*H0*BT1int;
%LMat is unchanged

%% Define operators in coupled spin basis
%SpinEx is already in coupled spin basis
Hdd2=BT2int*HddInt*BT2int';
E2=BT2int*E*BT2int';
H02=BT2int*H0Int*BT2int';
%LMat is unchanged

%% Define initial operators
SpinExPot=SpinEx;
HddPot=Hdd2;
EPot=E2;
H0Pot=H02;

Y=sqrt(PotentialFunc(rMin,Scale,LMat,SpinExPot,HddPot,H0Pot)-EPot).*IR; %Initial condition
if nargout>1
    unew = 1e-20*ones(NSubChannels,1);
    uout = unew;
    rout = rMin;
end

%% Calculate adiabatic potentials
r=linspace(rMin,rMax,5e3);
f=PotentialFunc(r,Scale,LMat,SpinExPot,HddPot,H0Pot)-repmat(EPot,[1,1,numel(r)]);
Adiabat=zeros(NSubChannels,numel(r));
for kk=1:numel(r)
    Adiabat(:,kk)=sort(eig(f(:,:,kk)));
end

%% Main Loop
BasisChangeFlag=false;
NumSegments=floor(rMax/BlockSize);
for jj=1:NumSegments
    %% Define integration segments due to memory limitations

    start_pos=rMin+(jj-1)*BlockSize;
    end_pos=start_pos+BlockSize;

    kmax=Adiabat(1,r>=start_pos & r<=end_pos);
    kmax=kmax(kmax<0);
    kmax=max(sqrt(abs(kmax)));
    if ~isempty(kmax)
        dr=IntParams.dr_scale*2*pi./kmax;
        dr=max(min(dr,IntParams.drMax),IntParams.drMin);
        dr=min(dr,BlockSize/2);
    end

    N=max(ceil((end_pos-start_pos)/dr),2);
    r2=linspace(start_pos,end_pos,N);
    %if N==1,error('Number of steps = 1');end;
    dr=diff(r2(1:2));
    h=dr/2;
    
    if nargout>1
        u = zeros(NSubChannels,N);
        u(:,1) = unew;
    end

    M=PotentialFunc(r2,Scale,LMat,SpinExPot,HddPot,H0Pot)-repmat(EPot,[1,1,N]);
    M2=PotentialFunc(r2+h,Scale,LMat,SpinExPot,HddPot,H0Pot)-repmat(EPot,[1,1,N]);

    %% Perform integration using Manolopoulos' log derivative method
    for nn=1:N-1
        p2=diag(M2(:,:,nn));
        p=sqrt(abs(p2));
        y1=p.*coth(p*h).*(p2>0)+p.*cot(p*h).*(p2<0)+1./h.*(p2==0);
        y2=p.*csch(p*h).*(p2>0)+p.*csc(p*h).*(p2<0)+1./h.*(p2==0);

        y1=diag(y1);
        y2=diag(y2);
        Mref=diag(p2);

        Qa=h/3*(M(:,:,nn)-Mref);
        Qc=4./h*((IR-h.^2/6.*(M2(:,:,nn)-Mref))\IR)-4./h*IR;
        Qb=h/3*(M(:,:,nn+1)-Mref);
        
        if nargout>1
            tmp = (y1+Qc)-(y2/(Y+y1+Qa))*y2;
            unew = (y2\(Y+y1+Qa))*unew;
            Y = (y1+Qb)-(y2/(tmp+y1+Qc))*y2;
            unew = (y2\(tmp+y1+Qc))*unew;
            if any(abs(unew)>1e30)
                u = u/norm(unew);
                uout = uout/norm(unew);
                unew = unew/norm(unew);
            end
            u(:,nn+1) = unew;
        else
            tmp=(y1+Qc)-(y2/(Y+y1+Qa))*y2;
            Y=(y1+Qb)-(y2/(tmp+y1+Qc))*y2;  
        end
    end
    
    if nargout>1
        uout = [uout u(:,2:end)]; %#ok<AGROW>
        rout = [rout r2(2:end)]; %#ok<AGROW>
    end

    if BasisChangeFlag==false && end_pos>=10
        Y=BT2int'*Y*BT2int;
        SpinExPot=SpinExInt;
        HddPot=HddInt;
        H0Pot=H0Int;
        EPot=E;
        BasisChangeFlag=true;
        if nargout>1
            unew = BT2int'*unew;
            for kk=1:size(uout,2)
                uout(:,kk) = BT2int'*uout(:,kk);
            end
        end
    end

end
b=r2(N);
% S-matrix MUST be calculated in basis where Hint (or k) is diagonal.
E=diag(E);
idx=E>=0;
E=diag(E);
k=sqrt(E(idx,idx));
LMatR2=LMat(idx,idx);
Y=Y(idx,idx);
Stmp=(k*RiccatiBessel(b,k,LMatR2,2,1)-RiccatiBessel(b,k,LMatR2,2,0)*Y)/(k*RiccatiBessel(b,k,LMatR2,1,1)-RiccatiBessel(b,k,LMatR2,1,0)*Y);
S=IR;
S(idx,idx)=Stmp;

if nargout>1
%     rout(end+1) = b;
    wf.r = rout;
    wf.u = uout;
end


end



function [Hint,BT1int,BT2int]=Generate_Hint(Hint0,BV1,BVint,BT21,B,gS_1,gS_2,gI_1,gI_2,Ahfs1,Ahfs2,NSpin1,NSpin2,mu_B,Scale,InitStateInt)

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
Etarget=Hint(InitStateInt,InitStateInt);    %in 1/Angstroms^2
Hint=Hint-Etarget*IR;
        
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
