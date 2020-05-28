classdef atomproperties < handle
    %ATOMPROPERTIES An abstract class for describing atomic ground state
    %properties
    properties(Abstract,Constant)
        mass    %Mass of the atom
        nspin   %Nuclear spin
        Ahfs    %Hyperfine constant in MHz
        gI      %Nuclear g-factor
        gS      %Electronic g-factor
    end
    
    methods
        function [E,V,Hint,label] = hyperfineSolve(self,B)
            %HYPERFINESOLVE Solves the ground-state hyperfine+zeeman
            %Hamiltonian
            %
            %   [E,V,Hint,label] = hyperfineSolve(self,B) solves the
            %   ground-state hyperfine+zeeman Hamiltonian for a particular
            %   magnetic field in Tesla.  E is a vector of eigenenergies, V
            %   is the transformation matrix from the eigenvectors to the
            %   uncoupled m_s and m_i states.  Hint is the total
            %   hyperfine+zeeman Hamiltonian.  label is the state labels in
            %   the uncoupled m_s, m_i basis
            S = 0.5;
            I = self.nspin;
            NI = 2*I+1;
            NS = 2*S+1;
            N = NI*NS;
            mS = -S:S;
            mI = -I:I;
            label = [reshape(repmat(mS,NI,1),N,1) repmat(mI(:),NS,1)];
            
            F = abs(I-S):abs(I+S);
            labelF = [];
            for nn=1:numel(F)
                mF = (-F(nn):F(nn))';
                labelF = [labelF;[repmat(F(nn),numel(mF),1) mF]]; %#ok<*AGROW>
            end
            
            %% Bare Hamilton
            H0 = zeros(N);
            for a=1:N
                mI1 = label(a,2);
                mS1 = label(a,1);
                for b=1:N
                    mI2 = label(b,2);
                    mS2 = label(b,1);
                    for c=1:numel(F)
                        if ~(abs(mI1+mS1)>F(c) || abs(mI2+mS2)>F(c) || (mI1+mS1)~=(mI2+mS2))
                            CG = ClebschGordan(I,S,F(c),mI1,mS1,mI1+mS1).*ClebschGordan(I,S,F(c),mI2,mS2,mI2+mS2);
                            K = F(c)*(F(c)+1)-I*(I+1)-S*(S+1);
                            tmp = CG*self.Ahfs/2*K;
                            H0(a,b) = H0(a,b)+tmp;
                        end
                    end
                end
            end
            
            %% Zeeman field
            HZ = diag(const.muB*B*(self.gS*label(:,1)+self.gI*label(:,2)))/const.h/1e6; %In [MHz]
            Hint = H0+HZ;
            [V,E] = eig(Hint,'vector');
            [E,idx] = sort(E);
            V = V(:,idx);
            E = diag(E);
        end
    end
    
    
end
    
 