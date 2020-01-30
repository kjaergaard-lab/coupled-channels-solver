classdef ScatteringMatrixElement < handle
    properties
        S
        T
        k
    end
    
    
    methods

        function obj = ScatteringMatrixElement(S,T)
            if nargin == 0
                obj.S = [];
                obj.T = [];
                return;
            else
                obj.S = S;
                obj.T = T;
            end
        end %end constructor()
        
        function r = a(obj)
            r = 1./(1i*obj.k).*(1-obj.S)./(1+obj.S);
        end
        
        
%         function B = subsref(obj,S)
%             if strcmp(S(1).type,'.') && strcmp(S(1).subs,'T')
%                 T = obj.S - obj.deltaValue;
%                 if numel(S) > 1
%                     B = subsref(T,S(2));
%                 else
%                     B = T;
%                 end
%             else
%                 B = builtin('subsref',obj,S);
%             end
%             
%         end
        
    end
    
end