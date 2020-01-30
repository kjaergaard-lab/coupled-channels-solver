classdef ScatteringMatrix < handle
    properties
        mass
        E
        B
        k
        
        BV
        symmetry
        SS
        targetIndex
    end
    
    methods
        function obj=ScatteringMatrix(Sin,BVin,symIn,InitStateLabel)
            if symIn == 0
                obj.BV = BVin;
                obj.symmetry = symIn;
                obj.targetIndex = find(all(obj.BV==repmat(InitStateLabel,size(Sin,1),1),2));
            elseif symIn == 1 || symIn == -1
                [Sin,obj.BV]=SymmetrizeSMatrixArbL(Sin,BVin,symIn);
                obj.symmetry = symIn;
                obj.targetIndex = find(all(obj.BV==repmat([InitStateLabel(1:2),sort(InitStateLabel(3:4),2)],size(Sin,1),1),2));
            else
                error('Input argument symIn must be either -1, 0, or 1');
            end
            
            obj.SS = ScatteringMatrixElement();
            S2 = shiftdim(Sin,2);
            for nn=1:size(Sin,1)
                S3 = squeeze(S2(:,:,obj.targetIndex));
                obj.SS(nn,1) = ScatteringMatrixElement(S3(:,nn),S3(:,nn)-(nn==obj.targetIndex));
            end
            
        end %end constructor()
        
        function set.E(obj,E)
            tmp = E(:)*const.kb/1e6;
            obj.E = E(:);
            obj.k = sqrt(2*obj.mass/const.hbar^2*tmp(:));  %#ok
            for nn=1:numel(obj.SS)
                obj.SS(nn).k = obj.k;
            end
        end
        
        function B = subsref(obj,S)
            switch S(1).type
                case '.'
                    switch S(1).subs
                        case 'S'
                            B = obj.SS(obj.targetIndex).S;
                        case 'T'
                            B = obj.SS(obj.targetIndex).T;
                        case 'a'
                            B = obj.SS(obj.targetIndex).a;
                        otherwise
                            B = builtin('subsref',obj,S);
                    end
                    
                case '()'
                    if numel(obj) > 1
                        B = obj(S(1).subs{1});
                    elseif numel(S(1).subs) == 1
                        B = obj.SS(S(1).subs{1});
                    elseif numel(S(1).subs) == 4
                        nn = find(all(obj.BV == repmat(cell2mat(S(1).subs),size(obj.BV,1),1),2));
                        if isempty(nn)
                            error('That state is not in the scattering matrix');
                        end
                        B = obj.SS(nn);
                    else
                        error('Not supported!');
                    end
                    
                    if numel(S) > 1
                        B = subsref(B,S(2:end));
                    end
                    
                otherwise
                    error('Index type not supported!');
            end
            
        end %end subsref()
        
        
        
    end
    
    
end %end class Scattering Matrix