function MultiChannel_Driver

% First run BasisSetGeneration('RbRb',0:2);

% E=[1:1:10,20:10:100,200:100:1000];
% B=linspace(645,675,2e2);

% E = [1e-3,1e-2];
% B = 300;

% B=930;
% E=linspace(1,1e3,1e2);

E = 300;
% B = unique([10e-3:10e-3:1,1:0.1:20,20]);
% B = linspace(913,914,3e2);
B = linspace(900,960,2e2);

DipoleFlag=1;   %Controls if dipole-dipole interactions are included.  Including it uses LOTS of memory
BasisSetFile='RbRb Basis Set L=0-2'; %Basis set - may need to run before
OutputFile='RbRb test';  %Output file
InitStateLabel=[0 0 1 1];         %In |L mL Int1 Int2> basis. States are in order of increasing energy

IntParams.rMin=1e-1;    %Starting distance in angstroms
IntParams.rMax=500;     %Final distance in angstroms
IntParams.dr_scale=5e-2;    %Step size scale factor -> approx linear relation with calculation time
IntParams.drMax=500;    %Maximum step size
IntParams.drMin=1e-3;   %Minimum step size
IntParams.ParSet=1;     %Use parallel processing

% for nn=1:numel(E)
%     nn
%     sm(nn)=MultiChannel(InitStateLabel,E(nn),B,OutputFile,BasisSetFile,DipoleFlag,IntParams);   %Note the empty string means output data is in variable Output only
% end

MultiChannel(InitStateLabel,E,B,OutputFile,BasisSetFile,DipoleFlag,IntParams);

% save(OutputFile);



