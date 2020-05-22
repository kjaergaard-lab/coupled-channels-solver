classdef scattoptions < handle
    %SCATTOPTIONS Options for solving coupled channels equations for
    %scattering states
    
    properties
        rmin        %Minimum integration distance in Angstroms
        rmax        %Maximum integration distance in Angstroms
        drscale     %Scale factor for calculating dr
        drmax       %Maximum integration step
        drmin       %Minimum integration step
        parallel    %Set to true (1) to use a parfor loop, false to use a normal for loop
        getwf       %Set to true to calculate wavefunctions
        dipole      %Set to use to include dipole-dipole interaction
        blocksize   %Size of integration blocks
    end
    
    methods
        function self = scattoptions(varargin)
            %SCATTOPTIONS Constructs a SCATTOPTIONS object with the given
            %properties
            %
            %   See also SCATTOPTIONS.SET
            self.setDefaults();
            self.set(varargin{:});
            if self.dipole
                self.blocksize = 5;
            end
        end
        
        function self = setDefaults(self)
            %SETDEFAULTS Sets the default values for this object
            self.rmin = 0.1;
            self.rmax = 500;
            self.drscale = 5e-2;
            self.drmax = 500;
            self.drmin = 1e-3;
            self.parallel = true;
            self.getwf = false;
            self.dipole = false;
            self.blocksize = 15;
        end
        
        function self = set(self,varargin)
            %SET Sets properties using name/value pairs.
            if mod(numel(varargin),2)~=0
                error('Must supply name/value pairs');
            else 
                for nn=1:2:numel(varargin)
                    if ~ischar(varargin{nn})
                        error('Optional argument %d must be a string',nn);
                    else
                        k = varargin{nn};
                        v = varargin{nn+1};
                        p = fieldnames(self);
                        for mm=1:numel(p)
                            if strcmpi(k,p{mm})
                                self.(p{mm}) = v;
                                break;
                            end
                        end
                    end
                end
            end
        end
        
        
    end
    
    
    
end