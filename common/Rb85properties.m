classdef Rb85properties < atomproperties
    %RB85PROPERTIES Defines properties related to Rb85 atoms
    properties(Constant)
        mass = scattconst.amu*84.911789738;  %[kg]
        nspin = 2.5;
        Ahfs = 3035.732439/(5/2+1/2);   %[MHz]
        gI = -0.0002936400;             
        gS = 2.00233113;
    end
end