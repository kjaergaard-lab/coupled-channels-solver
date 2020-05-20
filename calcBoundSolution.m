function [match,nodes,debugOut] = calcBoundSolution(r,Vfunc,E,ops,opt)

%% Parse arguments
if nargin<5
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

opt.direction = 1;
% opt.stopAtRoot = true;
opt.stopAfterR = true;
% opt.stopAtR = true;
% opt.stopR = 6;

optR = opt;
optR.stopAtRoot = false;
optR.stopAtR = true;
optR.stopAfterR = false;
optR.direction = -1;

if nargout>1 || opt.output
    opt.output = true;
    [ysL,rs,nL,rL,yL,zL,dbgL] = manolopoulos(r,Vfunc,E,ops,opt);
    optR.stopR = rs;
    optR.output = true;
    [ysR,~,nR,rR,yR,zR,dbgR] = manolopoulos(r,Vfunc,E,ops,optR);
else
    [ysL,rs,nL] = manolopoulos(r,Vfunc,E,ops,opt);
    optR.stopR = rs;
    [ysR,~,nR] = manolopoulos(r,Vfunc,E,ops,optR);
end


Ydiff = ysL-ysR;
[vecY,eigY] = eig(Ydiff,'vector');
[~,minIndex] = min(abs(eigY));
match = real(eigY(minIndex));
nodes = nR+nL;
wf = vecY(:,minIndex);

if nargout>2
%     u = [uL;flip(uR(1:end-1)/uR(end)*uL(end))];
    r = [rL;flip(rR(1:end-1))];
    Y = yL;
    Y(:,:,end+(1:(size(yR,3)-1))) = flip(yR(:,:,1:end-1),3);
%     debugOut.u = u;
    debugOut.r = r;
    debugOut.Y = Y;
    debugOut.yL = yL;
    debugOut.rL = rL;
    debugOut.eL = dbgL.eigOut;
    debugOut.yR = yR;
    debugOut.rR = rR;
    debugOut.eR = dbgR.eigOut;
    debugOut.wf = wf;
    debugOut.zL = zL;
    debugOut.zR = zR;
    debugOut.eigY = eigY;
    debugOut.Ydiff = Ydiff;
end

end