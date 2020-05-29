classdef K41properties < atomproperties
    %K41PROPERTIES Defines properties related to K41 atoms
    properties(Constant)
        mass = const.amu*40.96182576;  %[kg]
        nspin = 1.5;
        Ahfs = 254.013872/(3/2+1/2);   %[MHz]
        gI = -0.00007790600;             
        gS = 2.00229421;
    end
end