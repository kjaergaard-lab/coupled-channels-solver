function [match,nodes,debugOut] = calcBoundSolution(r,Vfunc,E,options)
narg = 4;

if nargin<narg || ~isfield(options,'stopAtRoot')
    options.stopAtRoot = false;
end

if nargin<narg || ~isfield(options,'stopR')
    options.stopR = 3;
end

if nargin<narg || ~isfield(options,'stopAfterR')
    options.stopAfterR = true;
end

options.direction = -1;

if nargout>1 || (isfield(options,'makeOutput') && options.makeOutput)
    options.makeOutput = true;
    [ysR,rs,rR,yR,uR] = manolopoulos(r,Vfunc,E,options);
    options.stopAtRoot = false;
    options.stopAtR = true;
    options.stopAfterR = false;
    options.stopR = rs;
    options.direction = 1;
    [ysL,~,rL,yL,uL] = manolopoulos(r,Vfunc,E,options);
else
    [ysR,rs] = manolopoulos(r,Vfunc,E,options);
    options.stopAtRoot = false;
    options.stopAtR = true;
    options.stopAfterR = false;
    options.stopR = rs;
    options.direction = 1;
    ysL = manolopoulos(r,Vfunc,E,options);
end

nodes = 0;
match = ysL-ysR;

if nargout>2
    u = [uL;flip(uR(1:end-1)/uR(end)*uL(end))];
    r = [rL;flip(rR)];
    Y = [yL;flip(yR)];
    debugOut.u = u;
    debugOut.r = r;
    debugOut.Y = Y;
    debugOut.uL = uL;
    debugOut.rL = rL;
    debugOut.uR = flip(uR/uR(end)*uL(end));
    debugOut.rR = flip(rR);
end

end