function [rout,blocks] = makegrid(E,ops,opt)
% MAKEGRID Creates a grid of inter-nuclear separations based on the energy,
% the potential, and the desired options
%
%   [rout,blocks] = MAKEGRID(E,ops,opt) Creates a grid rout separated into
%   blocks with constant step-sizes.  E is the energy to use, ops is an
%   instance of ALKALIOPERATORS, and opt is an instance of BOUNDOPTIONS
%
%   Relevant BOUNDOPTIONS properties are:
%
%   rmin:       minimum integration distance
%   rmax:       maximum integration distance
%   drscale:    scaling of the grid step relative to the shortest local wavelength
%   drmin:      minimum value of the step size
%   drmax:      maximum value of the step size
%   blocksize:  size of constant step-size blocks in Angstroms
if nargin<2
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

%% Create grid
numSegments = ceil(opt.rmax/opt.blocksize);
rtmp = linspace(opt.rmin,opt.rmin+numSegments*opt.blocksize,1e3);
f = ops.potential(rtmp,E);
adiabat = zeros(ops.Nch,numel(rtmp));
for kk=1:numel(rtmp)
    adiabat(:,kk) = sort(eig(f(:,:,kk)));
end

rout = [];
blocks = zeros(numSegments,2);

for jj=1:numSegments
    rStart = opt.rmin+(jj-1)*opt.blocksize;
    rEnd = rStart+opt.blocksize;
    kmax = adiabat(1,rtmp>=rStart & rtmp<=rEnd);
    kmax = max(sqrt(abs(kmax(kmax<0))));
%     kmax = max(sqrt(abs(kmax)));
    if ~isempty(kmax)
        dr = opt.drscale*2*pi/kmax;
        dr = max(min(dr,opt.drmax),opt.drmin);
        dr = min(dr,opt.blocksize/2);
    end
    
    N = max(ceil(abs(rEnd-rStart)/dr),2);
    r = linspace(rStart,rEnd,N)';
    blocks(jj,:) = numel(rout)+[1,N];
    rout = [rout;r(1:end-1)]; %#ok<AGROW>
    
end

rout(end+1) = r(end);

end