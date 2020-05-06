function [rout,Vout] = makegrid(Vfunc,E,opt)

if nargin<3
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

%% Create grid
numSegments = ceil(opt.rmax/opt.blocksize);
rtmp = linspace(opt.rmin,numSegments*opt.blocksize,5e3);
adiabat = Vfunc(rtmp)-E;

rout = [];

for jj=1:numSegments
    rStart = opt.rmin+(jj-1)*opt.blocksize;
    rEnd = rStart+opt.blocksize;
    kmax = adiabat(rtmp>=rStart & rtmp<=rEnd);
    kmax = max(sqrt(abs(kmax)));
    if ~isempty(kmax)
        dr = opt.drscale*2*pi/kmax;
        dr = max(min(dr,opt.drmax),opt.drmin);
        dr = min(dr,opt.blocksize/2);
    end
    
    N = max(ceil(abs(rEnd-rStart)/dr),2);
    r = linspace(rStart,rEnd,N)';
    rout = [rout;r(1:end-1)]; %#ok<AGROW>
    
end

rout(end+1) = r(end);

Vout = Vfunc(rout);

end