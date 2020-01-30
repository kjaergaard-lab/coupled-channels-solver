function V=RbRbPotentialDD(x,Scale,LMat,SpinEx,Hdd,H0)
%Calculates r-dependent potentials.  
a_Bohr=0.52917721092;               %Angstroms

%% Dipole-dipole terms
aSO=-0.0416./a_Bohr.^3;
R_SO=7.5*a_Bohr;
b=0.7196./a_Bohr;

%% Hyperfine constant variation
cf=-0.0778;
R0=11*a_Bohr;
DeltaR=0.5*a_Bohr;

Nx=numel(x);
N=size(LMat,1);
x=reshape(x(:),[1,1,Nx]);
LMat=repmat(LMat.*(LMat+eye(N)),[1,1,Nx]);  %Matrix of L values

%V=VDirect + SpinEx * J
VDirect=0.25*(3*RbBOPotentialsv2(x,1)+RbBOPotentialsv2(x,0)); %Direct term
J=(RbBOPotentialsv2(x,1)-RbBOPotentialsv2(x,0));  %Exchange term
VDirect=repmat(VDirect,[N,N,1]).*repmat(eye(N),[1,1,Nx]);
J=repmat(J,[N,N,1]).*repmat(SpinEx,[1,1,Nx]);

%Dipole-dipole term
Hdd=Scale*bsxfun(@times,Hdd,1./x.^3+aSO.*exp(-b*(x-R_SO)));

%Variation of hyperfine constants
H0=bsxfun(@times,H0,cf./(exp((x-R0)./DeltaR)+1));

%Total potential
V=Scale*(VDirect+J)+Hdd+bsxfun(@times,LMat,x.^(-2))+H0;
end