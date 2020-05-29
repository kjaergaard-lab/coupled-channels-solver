classdef alkalioperators
    %ALKALIOPERATORS Operators for use with coupled-channels solvers
    %
    %   Contains operators used in calculating potential operators for
    %   solving coupled-channels equations for alkali metal atoms
    %
    properties
        Nch         %Number of channels
        basis       %0 if internal, 1 if BV2
        pair        %Atom pair string indicating which potential to use
        scale       %Scale factor
        L           %Matrix of L values
        SpinProj    %Matrix for spin projections
        Hdd         %Dipole-dipole operator
        H0          %Internal spin-nucleus Hamiltonian
        Hint        %Internal Hamiltonian with Zeeman Hamiltonian
        refE        %Reference energy level
        
        BT2int      %Transformation matrix from internal states to BV2
        U           %Transformation matrix
        
    end
    
    methods
        function self = alkalioperators(pair,scale,L,SpinProj,Hdd,H0,Hint,refE,BT2int)
            %ALKALIOPERATORS Constructs an instance of the class
            %
            %   self = ALKALIOPERATORS(pair,scale,L,SpinProj,Hdd,H0,Hint,refE,BT2int)
            %   constructs an instance of the class with the given
            %   parameters.  Pair is a character string of either 'rbrb',
            %   'kk', or 'krb' which denotes which potential function to
            %   use.  'scale' is the conversion from cm^{-1} to inverse
            %   wavenumbers in Angstroms^{-2}.  L is the matrix of L
            %   values.  SpinProj is the matrix representation of S1\cdot S2
            %   Hdd is the dipole-dipole operator.  H0 is the bare hyperfine
            %   Hamiltonian without the Zeeman effect.  Hint is the total
            %   internal Hamiltonian with the Zeeman effect.  refE is the
            %   reference energy - i.e. it is the asymptotic energy of the
            %   entrance channel.  BT2int is the transformation matrix
            %   taking vectors in the internal basis to the electron spin
            %   coupled basis.
            %
            %   All operators MUST be in basis 2 (electron spin coupled)
            self.pair = pair;
            self.scale = scale;
            self.L = L;
            self.SpinProj = SpinProj;
            self.Hdd = Hdd;
            self.H0 = H0;
            self.Hint = Hint;
            self.BT2int = BT2int;
            self.Nch = size(self.L,1);
            self.basis = 1;
            self.U = self.BT2int;
            self.refE = refE;
        end
        
        function self = rotate(self)
            %ROTATE Rotates the basis from spin-coupled to internal or vice
            %versa depending on which basis the operators are currently in
            if self.basis == 0
                self = self.int2spin;
            elseif self.basis == 1
                self = self.spin2int;
            end
        end
        
        function self = int2spin(self)
            %INT2SPIN Changes the basis from the internal basis to the
            %spin-coupled basis
            self.U = self.BT2int;
            if self.basis == 1
                return;
            end
            self.SpinProj = self.U*self.SpinProj*self.U';
            self.Hdd = self.U*self.Hdd*self.U';
            self.H0 = self.U*self.H0*self.U';
            self.Hint = self.U*self.Hint*self.U';
            self.basis = 1;
        end
        
        function self = spin2int(self)
            %SPIN2INT Changes the basis from the spin-coupled basis to the
            %internal basis
            self.U = self.BT2int';
            if self.basis == 0
                return;
            end
            self.SpinProj = self.U*self.SpinProj*self.U';
            self.Hdd = self.U*self.Hdd*self.U';
            self.H0 = self.U*self.H0*self.U';
            self.Hint = self.U*self.Hint*self.U';
            self.basis = 0;
        end
        
        function V = potential(self,r,E)
            %POTENTIAL Calculates the potential operator for a given vector
            %of positions r and an energy E above the reference energy
            if strcmpi(self.pair,'krb')
                V = KRbPotentialDD(r,self.scale,self.L,self.SpinProj,self.Hdd,self.H0)+self.Hint-(E+self.refE)*eye(self.Nch);
            elseif strcmpi(self.pair,'kk')
                V = KKPotentialDD(r,self.scale,self.L,self.SpinProj,self.Hdd,self.H0)+self.Hint-(E+self.refE)*eye(self.Nch);
            elseif strcmpi(self.pair,'rbrb')
                V = RbRbPotentialDD(r,self.scale,self.L,self.SpinProj,self.Hdd,self.H0)+self.Hint-(E+self.refE)*eye(self.Nch);
            end
        end
        
    end
    
    
    
    
end