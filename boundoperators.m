classdef boundoperators
    properties
        Nch         %Number of channels
        basis       %0 if internal, 1 if BV2
        scale       %Scale factor
        L           %Matrix of L values
        SpinProj    %Matrix for spin projections
        Hdd         %Dipole-dipole operator
        H0          %Internal spin-nucleus Hamiltonian
        Hint        %Internal Hamiltonian with Zeeman Hamiltonian
        
        BT2int      %Transformation matrix from internal states to BV2
        U           %Transformation matrix
        
    end
    
    methods
        function op = boundoperators(scale,L,SpinProj,Hdd,H0,Hint,BT2int)
            op.scale = scale;
            op.L = L;
            op.SpinProj = SpinProj;
            op.Hdd = Hdd;
            op.H0 = H0;
            op.Hint = BT2int*Hint*BT2int';
            op.BT2int = BT2int;
            op.Nch = size(op.L,1);
            op.basis = 1;
            op.U = op.BT2int;
        end
        
        function op = rotate(op)
            if op.basis == 0
                op = op.int2spin;
            elseif op.basis == 1
                op = op.spin2int;
            end
        end
        
        function op = int2spin(op)
            op.U = op.BT2int;
            if op.basis == 1
                return;
            end
            op.SpinProj = op.U*op.SpinProj*op.U';
            op.Hdd = op.U*op.Hdd*op.U';
            op.H0 = op.U*op.H0*op.U';
            op.Hint = op.U*op.Hint*op.U';
            op.basis = 1;
        end
        
        function op = spin2int(op)
            op.U = op.BT2int';
            if op.basis == 0
                return;
            end
            op.SpinProj = op.U*op.SpinProj*op.U';
            op.Hdd = op.U*op.Hdd*op.U';
            op.H0 = op.U*op.H0*op.U';
            op.Hint = op.U*op.Hint*op.U';
            op.basis = 0;
        end
        
    end
    
    
    
    
end