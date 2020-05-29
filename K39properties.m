classdef K39properties < atomproperties
    %K39PROPERTIES Defines properties related to K39 atoms
    properties(Constant)
        mass = const.amu*38.96370668;  %[kg]
        nspin = 1.5;
        Ahfs = 461.7197202/(3/2+1/2);  %[MHz]
        gI = -0.00014193489;             
        gS = 2.00229421;
    end
end