classdef atomproperties < handle
    properties(Abstract,Constant)
        mass
        nspin
        Ahfs
        gI
        gS
    end
    
    methods
        function [E,V,Hint,label] = hyperfineSolve(self,B)
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
    
 