classdef boundresults < handle
    properties
        Nbound
        basis
        E
        r
        wf
        BV2
        BVint
        BT2int
    end
    
    methods
        function self = boundresults(E,r,wf)
            if nargin>0
                self.E = E;
                self.r = r;
                self.wf = wf;
                self.Nbound = numel(E);
            end
        end
        
        function N = count(self)
            N = numel(self.E);
            self.Nbound = N;
        end
        
        function rotate(self)
            if self.basis == 0
                U = self.BT2int;
                self.basis = 1;
            else
                U = self.BT2int';
                self.basis = 0;
            end
            for nn=1:numel(self.wf)
                for mm=1:size(self.wf{nn},2)
                    self.wf{nn}(:,mm) = U*self.wf{nn}(:,mm);
                end
            end
                
        end
        
        
    end
    
end