classdef alkalioperators
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
            if self.basis == 0
                self = self.int2spin;
            elseif self.basis == 1
                self = self.spin2int;
            end
        end
        
        function self = int2spin(self)
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