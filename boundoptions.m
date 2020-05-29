classdef boundoptions
    properties
        stopAtRoot  %Stop integration at a root in the log-derivative
        stopAtR     %Stop at a particular value of r
        stopAfterR  %Stop integration on condition after a value of r
        stopR       %Value of r to use in stop conditions
        
        direction   %Direction to integrate (positive/negative values)
        output      %Create output vectors
        
        iter        %Number of iterations of the bound-state energy finder
        rangeiter   %Number of iterations for finding energy ranges
        tolF        %Tolerance for the bound-state energy finder
        tolE        %Tolerance in energy for the bound-state energy finder
        debug       %Create debugging outputs
        pauseDelay  %Delay during debugging
        
        rmin        %Minimum r to use
        rmax        %Maximum r to use
        drmin       %Minimum dr
        drmax       %Maximum dr
        drscale     %Scaling to use when calculating dr
        blocksize   %Size of blocks with which to calculate dr
        blocks      %Indices of blocks
        
        usedipole   %Include dipole-dipole interaction
        changeR     %R value at which to change bases
        
    end
    
    methods
        function opt = boundoptions(varargin)
            opt = opt.setDefaults;
            opt = opt.set(varargin{:});
        end
        
        function opt = setDefaults(opt)
            opt.stopAtRoot = false;
            opt.stopAtR = false;
            opt.stopAfterR = false;
            opt.stopR = 6;
            
            opt.direction = 1;
            opt.output = false;
            
            opt.iter = 10;
            opt.rangeiter = 4;
            opt.tolF = 1e-6;
            opt.tolE = 1e-6;
            opt.debug = false;
            opt.pauseDelay = 0.1;
            
            opt.rmin = 2;
            opt.rmax = 500;
            opt.drmin = 1e-3;
            opt.drmax = 50;
            opt.drscale = 1e-2;
            opt.blocksize = 15;
            
            opt.usedipole = false;
            opt.changeR = 10;
        end
        
        function opt = set(opt,varargin)
            if mod(numel(varargin),2)~=0
                error('Must supply name/value pairs');
            else 
                for nn=1:2:numel(varargin)
                    if ~ischar(varargin{nn})
                        error('Optional argument %d must be a string',nn);
                    else
                        k = varargin{nn};
                        v = varargin{nn+1};
                        p = fieldnames(opt);
                        for mm=1:numel(p)
                            if strcmpi(k,p{mm})
                                opt.(p{mm}) = v;
                                break;
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    
    
end