%Generates constants
%% Universal constants
amu=1.660538921e-27;
% mK = 39.96399848*1.660538921e-27;   %kg
mK = 39.963998166*amu;  %Mass of 40K
mRb = 86.909180527*amu;    %Mass of 87Rb
mK41=40.96182576*amu;   %Mass of 41K
mRb85 = 84.911789738*amu;
mK39 = 38.96370668*amu;

% mRb = 1.44316060e-25;
hbar=6.62606957e-34/(2*pi);         %Js
kb=1.3806488e-23;                  %J/K
LightSpeed=299792458;               %m/s
mu_B=9.27400968e-24;                %J/T
a_Bohr=0.52917721092;               %Angstroms
mu0=4*pi*1e-7;                      %N/A^2
gS=2.0023193043622;
FineStructure=1/137.035999074;
EHartree=hbar*LightSpeed*FineStructure/(a_Bohr*1e-10);

%% Rubidium-87 hyperfine parameters
% Ahfs_Rb=6834.682610904/2;   %MHz
Ahfs_Rb=3417.34130545215;
Ahfs_Rb=Ahfs_Rb*1e6/LightSpeed*1e-2;    %cm^{-1}
gI_Rb=-0.0009951414;
gS_Rb=2.00233113;
% gS_Rb=gS;

%% Rubidium-85 hyperfine parameters
% Ahfs_Rb=6834.682610904/2;   %MHz
Ahfs_Rb85=3035.732439/(5/2+1/2);
Ahfs_Rb85=Ahfs_Rb85*1e6/LightSpeed*1e-2;    %cm^{-1}
gI_Rb85=-0.0002936400;
gS_Rb85=2.00233113;
% gS_Rb=gS;

%% Potassium-40 hyperfine parameters
Ahfs_K=-285.7308;           %MHz
Ahfs_K=Ahfs_K*1e6/LightSpeed*1e-2;    %cm^{-1}
gI_K=0.000176490;
gS_K=2.00229421;
% gS_K=gS;

%% Potassium-41 hyperfine parameters
Ahfs_K41=254.013872/(3/2+1/2);           %MHz
Ahfs_K41=Ahfs_K41*1e6/LightSpeed*1e-2;    %cm^{-1}
gI_K41=-0.00007790600;
gS_K41=2.00229421;
% gS_K=gS;

%% Potassium-39 hyperfine parameters
Ahfs_K39=461.7197202/(3/2+1/2);           %MHz
Ahfs_K39=Ahfs_K39*1e6/LightSpeed*1e-2;    %cm^{-1}
gI_K39=-0.00014193489;
gS_K39=2.00229421;
% gS_K=gS;