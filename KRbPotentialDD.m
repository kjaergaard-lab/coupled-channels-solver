function V=KRbPotentialDD(x,Scale,LMat,SpinEx,Hdd,H0)
%Calculates r-dependent potentials.  
a_Bohr=0.52917721092;               %Angstroms

%% Dipole-dipole terms
aSO=-0.013./a_Bohr.^3;
R_SO=7.5*a_Bohr;
b=0.7196./a_Bohr;

%% Hyperfine constant variation (none specified for KRb)
% cf=0;
% R0=0;
% DeltaR=0;

%% Calculate values
Nx=numel(x);
N=size(LMat,1);
x=reshape(x(:),[1,1,Nx]);
LMat=repmat(LMat.*(LMat+eye(N)),[1,1,Nx]);  %Matrix of L values

%V=VDirect + SpinEx * J
VDirect=0.25*(3*KRbBOPotentialsv3(x,1)+KRbBOPotentialsv3(x,0)); %Direct term
J=(KRbBOPotentialsv3(x,1)-KRbBOPotentialsv3(x,0));  %Exchange term
VDirect=repmat(VDirect,[N,N,1]).*repmat(eye(N),[1,1,Nx]);
J=repmat(J,[N,N,1]).*repmat(SpinEx,[1,1,Nx]);

%Dipole-dipole term
Hdd=bsxfun(@times,Hdd,1./x.^3+aSO.*exp(-b*(x-R_SO)));

%Variation of hyperfine constants
% H0=bsxfun(@times,H0,cf./(exp((x-R0)./DeltaR)+1));
H0=0;

%Total potential
V=Scale*(VDirect+J)+Hdd+bsxfun(@times,LMat,x.^(-2))+H0;
end