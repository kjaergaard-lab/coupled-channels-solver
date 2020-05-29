classdef Rb87properties < atomproperties
    %RB87PROPERTIES Defines properties related to Rb87 atoms
    properties(Constant)
        mass = const.amu*86.909180527;  %[kg]
        nspin = 1.5;
        Ahfs = 3417.34130545215;        %[MHz]
        gI = -0.0009951414;             
        gS = 2.00233113;
    end
end