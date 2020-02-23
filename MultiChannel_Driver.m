function MultiChannel_Driver

% First run BasisSetGeneration('RbRb',0:2);


E = [1e-3,1e-2,1e-1,1:20,25:5:100,100:10:1e3];
B = unique([linspace(911,913,2e2),linspace(913,918,1e2)]);

DipoleFlag=1;   %Controls if dipole-dipole interactions are included.  Including it uses LOTS of memory
BasisSetFile='RbRb Basis Set L=0-2'; %Basis set - may need to run before
%OutputFile='RbRb test';  %Output file
InitStateLabel=[2 0 1 1];         %In |L mL Int1 Int2> basis. States are in order of increasing energy

IntParams.rMin=1e-1;    %Starting distance in angstroms
IntParams.rMax=500;     %Final distance in angstroms
IntParams.dr_scale=5e-2;    %Step size scale factor -> approx linear relation with calculation time
IntParams.drMax=500;    %Maximum step size
IntParams.drMin=1e-3;   %Minimum step size
IntParams.ParSet=1;     %Use parallel processing

for nn=1:numel(E)
    formatSpec = './data/RbRb 1-1 d-wave at 930 G with DD, E = %d %s';
    if E(nn) < 1
       OutputFile = sprintf(formatSpec,E(nn)*1e3,'nK');
    else
       OutputFile = sprintf(formatSpec,E(nn),'uK');
    end
    fprintf(1,'Run: %d/%d, File: %s\n',nn,numel(E),OutputFile);
    MultiChannel(InitStateLabel,E(nn),B,OutputFile,BasisSetFile,DipoleFlag,IntParams);   %Note the empty string means output data is in variable Output only
end

%MultiChannel(InitStateLabel,E,B,OutputFile,BasisSetFile,DipoleFlag,IntParams);

% save(OutputFile);



