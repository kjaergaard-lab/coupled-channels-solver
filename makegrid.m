function [rout,Vout] = makegrid(Vfunc,E,varargin)
%% Default values
rmin = 0.1;
rmax = 500;
drmin = 1e-3;
drmax = 50;
drscale = 1e-2;
blocksize = 15;

%% Parse arguments
if numel(varargin)>0 && mod(numel(varargin),2)~=0
    error('Must supply name/value pairs for additional arguments');
elseif numel(varargin)==1
    options = flattenstruct(varargin{1});
else
    options = varargin;
end

for nn=1:2:numel(options)
    if ~ischar(options{nn})
        error('Optional argument %d is not a string',nn);
    else
        v = options{nn+1};
        switch lower(options{nn})
            case 'drmin'
                drmin = v;
            case 'drmax'
                drmax = v;
            case 'drscale'
                drscale = v;
            case 'rmin'
                rmin = v;
            case 'rmax'
                rmax = v;
            case 'blocksize'
                blocksize = v;
        end
    end
end

%% Create grid
numSegments = ceil(rmax/blocksize);
rtmp = linspace(rmin,numSegments*blocksize,5e3);
adiabat = Vfunc(rtmp)-E;

rout = [];

for jj=1:numSegments
    rStart = rmin+(jj-1)*blocksize;
    rEnd = rStart+blocksize;
    kmax = adiabat(rtmp>=rStart & rtmp<=rEnd);
    kmax = max(sqrt(abs(kmax)));
    if ~isempty(kmax)
        dr = drscale*2*pi/kmax;
        dr = max(min(dr,drmax),drmin);
        dr = min(dr,blocksize/2);
    end
    
    N = max(ceil(abs(rEnd-rStart)/dr),2);
    r = linspace(rStart,rEnd,N)';
    rout = [rout;r(1:end-1)]; %#ok<AGROW>
    
end

rout(end+1) = r(end);

Vout = Vfunc(rout);

end