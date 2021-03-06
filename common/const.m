classdef const < handle
    properties(Constant)
        %% Universal constants
        amu=1.660538921e-27;
        mElectron = 9.10938291e-31;
        mK = 39.963998166*const.amu;  %Mass of 40K
        mRb = 86.909180527*const.amu;    %Mass of 87Rb
        mKRb = const.mK*const.mRb/(const.mRb+const.mK);
        mK41=40.96182576*const.amu;   %Mass of 41K
        mRb85 = 84.911789738*const.amu;
        mK39 = 38.96370668*const.amu;

        e = 1.609e-19;    %C
        h = 6.62606957e-34;
        hbar = const.h/(2*pi);         %Js
        kb = 1.3806488e-23;      %J/K
        c = 299792458;           %m/s
        muB = 9.27400968e-24;   %J/T
        aBohr = 0.52917721092e-10;    %m
        mu0 = 4*pi*1e-7;         %N/A^2
        eps0 = 1./(const.mu0*const.c^2);
        gS = 2.0023193043622;
        alpha = 1/137.035999074;
        EHartree = const.hbar*const.c*const.alpha/const.aBohr;
        
        g = 9.81;   %m/s^2
        
        Rb87 = Rb87properties;
        Rb85 = Rb85properties;
        K39 = K39properties;
        K40 = K40properties;
        K41 = K41properties;
    end
    
    methods(Static)
        function a=cm2K
            %Converts an energy given in cm^{-1} to K
            a = 1e2*const.h*const.c/const.kb;
        end
        
        function a = K2A(mass)
            %K2A Converts an energy in Kelvin to a wavenumber^2 in
            %Angstrom^2.
            %   A = K2A(MASS) calculates conversion constant using given
            %   reduced mass
            a = 2*mass/const.hbar^2*1e-20*const.kb;
        end
        
    end
end
