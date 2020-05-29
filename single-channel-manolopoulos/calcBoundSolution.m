function [match,nodes,debugOut] = calcBoundSolution(r,Vfunc,E,opt)

%% Parse arguments
if nargin<4
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

opt.direction = 1;
opt.stopAtRoot = true;

opt2 = opt;
opt2.stopAtRoot = false;
opt2.stopAtR = true;
opt2.stopAfterR = false;
opt2.direction = -1;

if nargout>1 || opt.output
    opt.output = true;
    [ysL,rs,nL,rL,yL,uL] = manolopoulos(r,Vfunc,E,opt);
    opt2.stopR = rs;
    [ysR,~,nR,rR,yR,uR] = manolopoulos(r,Vfunc,E,opt2);
else
    [ysL,rs,nL] = manolopoulos(r,Vfunc,E,opt);
    opt2.stopR = rs;
    [ysR,~,nR] = manolopoulos(r,Vfunc,E,opt2);
end

nodes = nR+nL;
match = ysL-ysR;

if nargout>2
    u = [uL;flip(uR(1:end-1)/uR(end)*uL(end))];
    r = [rL;flip(rR(1:end-1))];
    Y = [yL;flip(yR(1:end-1))];
    debugOut.u = u;
    debugOut.r = r;
    debugOut.Y = Y;
    debugOut.uL = uL;
    debugOut.rL = rL;
    debugOut.uR = flip(uR/uR(end)*uL(end));
    debugOut.rR = flip(rR);
end

end