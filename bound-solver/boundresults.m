classdef boundresults < handle
% BOUDNRESULTS Provides a class for describing results from the bound state
% solver
    properties
        Nbound  %Number of bound states
        basis   %Basis of the wavefunctions
        E       %Bound state energies as an 1xNbound cell array
        r       %Vector of internuclear separations
        wf      %Wavefunctions as a 1xNbound cell array for each eigenvector
        bv2     %Basis vectors in the electron spin coupled basis
        bvint   %Basis vectors in the internal basis
        bt2int  %Basis transformation matrix from internal basis to spin coupled basis
    end
    
    methods
        function self = boundresults(E,r,wf)
            %BOUNDRESULTS Constructs a BOUNDRESULTS object
            %
            %   res = BOUNDRESULTS() constructs an object with no
            %   properties set
            %
            %   res = BOUNDRESULTS(E,r,wf) constructs an object with bound
            %   state energies E, separations r, and eigenvectors wf
            if nargin>0
                self.E = E;
                self.r = r;
                self.wf = wf;
                self.Nbound = numel(E);
            end
        end
        
        function N = count(self)
            %COUNT Counts the number of bound states
            %
            %   N = self.count() Sets Nbound property and returns the
            %   result
            N = numel(self.E);
            self.Nbound = N;
        end
        
        function rotate(self)
            %ROTATE Rotates the basis of the wavefunctions from one basis
            %to the other
            %
            %   self.rotate() rotates the wavefunctions from the basis they
            %   are in to the other basis using the objects value of bt2int
            if self.basis == 0
                U = self.bt2int;
                self.basis = 1;
            else
                U = self.bt2int';
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