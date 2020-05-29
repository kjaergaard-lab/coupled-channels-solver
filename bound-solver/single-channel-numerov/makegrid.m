function [x,V] = makegrid(xs,Vfunc,E,options)
if nargin<4 || ~isfield(options,'dxmin')
    options.dxmin = 1e-3;
end

if nargin<4 || ~isfield(options,'dxmax')
    options.dxmax = 100;
end

if nargin<4 || ~isfield(options,'dxscale')
    options.dxscale = 5e-2;
end

if nargin<4 || ~isfield(options,'initsize')
    options.initsize = 1e4;
end

if nargin<4 || ~isfield(options,'blocksize')
    options.blocksize = 1;
end


Vs = Vfunc(xs,E);
[zi,zdi] = findzeros(Vs);
if isempty(zi)
    error('No zeros found');
elseif numel(zi) == 1
    zi2 = propwkb(xs,Vs,zi,zdi);
    if zdi == -1
        xstart = xs(zi2);
        xend = xs(end);
    else
        xstart = xs(1);
        xend = xs(zi2);
    end
elseif numel(zi) == 2
    zi2 = zi;
    for nn=1:numel(zi)
        zi2(nn) = propwkb(xs,Vs,zi(nn),zdi(nn));
    end
    xstart = xs(min(zi2));
    xend = xs(max(zi2));
else
    error('More than two zeros found.  Unable to make grid');
end

dx = zeros(options.initsize,1);
x = xstart;
nn = 0;
while x<xend
    dxnew = options.dxscale/sqrt(abs(Vfunc(x,E)));
    if dxnew<options.dxmin
        dxnew = options.dxmin;
    elseif dxnew>options.dxmax
        dxnew = options.dxmax;
    end
    if nn == numel(dx)
        dx = [dx;zeros(options.initsize,1)]; %#ok<AGROW>
    end
    nn = nn+1;
    dx(nn) = dxnew;
    x = x+dxnew;
end
dx = dx(1:nn);

if options.blocksize == 1
    x = cumsum(dx);
elseif options.blocksize > 1
    numBlocks = ceil(numel(dx)/options.blocksize);
    for nn=1:numBlocks
        idx = ((nn-1)*options.blocksize+1):min(nn*options.blocksize,numel(dx));
        dx(idx) = min(dx(idx)); 
    end
    x = cumsum(dx);
else
    error('Block size must be >= 1');
end

x = x(x<=xend);
V = Vfunc(x,0);

end