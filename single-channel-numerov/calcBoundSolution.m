function [match,u,nodes,debugOut] = calcBoundSolution(x,W,options)

N = numel(x);
options.stopAtMax = 1;
options.stopXBelow = 3;
[uR,duR,~,idx,nR] = numerov(x,[1e-20,1.1e-20],W,-1,options);
[uL,duL,~,~,nL] = numerov(x(1:(N-idx+3)),[0,1e-20],W(1:(N-idx+3)),1);

nodes = nL+nR;
match = duL(end)./uL(end-1)-duR(end)./uR(end-1);
if nargout>1
    u = [uL(1:end-2);flip(uR(1:end-1)/uR(end-1)*uL(end-1))];
    u = u./sqrt(sum(u.^2*abs(x(2)-x(1))));
end

if nargout>3
    debugOut.uL = uL;
    debugOut.duL = duL;
    debugOut.xL = x(1:numel(uL));
    debugOut.uR = flip(uR/uR(end-1)*uL(end-1));
    debugOut.duR = flip(duR);
    debugOut.xR = x((end-idx+1):end);
end

end