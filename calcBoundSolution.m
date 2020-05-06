function [match,nodes,debugOut] = calcBoundSolution(r,Vfunc,E,opt)

%% Parse arguments
if nargin<4
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

opt.direction = -1;
opt = opt.set('stopAfterR',true);

if nargout>1 || opt.output
    opt.output = true;
    [ysR,rs,nR,rR,yR,uR] = manolopoulos(r,Vfunc,E,opt);
    opt.stopAtRoot = false;
    opt.stopAtR = true;
    opt.stopAfterR = false;
    opt.stopR = rs;
    opt.direction = 1;
    [ysL,~,nL,rL,yL,uL] = manolopoulos(r,Vfunc,E,opt);
else
    [ysR,rs,nR] = manolopoulos(r,Vfunc,E,opt);
    opt.stopAtRoot = false;
    opt.stopAtR = true;
    opt.stopAfterR = false;
    opt.stopR = rs;
    opt.direction = 1;
    [ysL,~,nL] = manolopoulos(r,Vfunc,E,opt);
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