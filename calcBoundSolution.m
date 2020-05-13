function [match,nodes,debugOut] = calcBoundSolution(r,Vfunc,E,opt)

%% Parse arguments
if nargin<4
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

opt.direction = 1;
opt.stopAtRoot = true;
% opt.stopR = 6;

optR = opt;
optR.stopAtRoot = false;
optR.stopAtR = true;
optR.stopAfterR = false;
optR.direction = -1;

if nargout>1 || opt.output
    opt.output = true;
    [ysL,rs,nL,rL,yL,zL] = manolopoulos(r,Vfunc,E,opt);
    optR.stopR = rs;
    [ysR,~,nR,rR,yR,zR] = manolopoulos(r,Vfunc,E,optR);
else
    [ysL,rs,nL] = manolopoulos(r,Vfunc,E,opt);
    optR.stopR = rs;
    [ysR,~,nR] = manolopoulos(r,Vfunc,E,optR);
end


Ydiff = ysL-ysR;
[vecY,eigY] = eig(Ydiff);
[~,minIndex] = min(abs(diag(eigY)));
match = eigY(minIndex,minIndex);
nodes = nR+nL+sum(diag(eigY)<0);
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
    debugOut.yR = yR;
    debugOut.rR = rR;
    debugOut.wf = wf;
    debugOut.zL = zL;
    debugOut.zR = zR;
    debugOut.eigY = eigY;
    debugOut.Ydiff = Ydiff;
end

end