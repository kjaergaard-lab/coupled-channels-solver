classdef K40properties < atomproperties
    %K40PROPERTIES Defines properties related to K40 atoms
    properties(Constant)
        mass = scattconst.amu*39.963998166;  %[kg]
        nspin = 4;
        Ahfs = -285.7308;               %[MHz]
        gI = 0.000176490;             
        gS = 2.00229421;
    end
end