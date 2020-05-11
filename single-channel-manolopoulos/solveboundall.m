function [Eout,u,E] = solvebound(x,V,Erange,options)

if nargin<4 || ~isfield(options,'maxIter')
    options.maxIter = 100;
end
if nargin<4 || ~isfield(options,'tol')
    options.tol = 1e-5;
end
if nargin<4 || ~isfield(options,'debug')
    options.debug = 0;
end
if nargin<4 || ~isfield(options,'pauseTime')
    options.pauseTime = 0.1;
end

err = 1;

x = x(:);
N = numel(x);
xflip = flip(x);
Vflip = flip(V);
options_numerov.stopAtMax = 1;
solCount = 0;

match = zeros(options.maxIter,1);
E = zeros(options.maxIter,1);
E(1:2) = Erange;
for nn=1:options.maxIter
    %% Generate candidate solution
    [uR,duR,~,idx,ncR] = numerov(xflip,[1e-20,1.1e-20],Vflip-E(nn),options_numerov);
    [uL,duL,~,~,ncL] = numerov(x(1:(N-idx+3)),[0,1e-20],V(1:(N-idx+3))-E(nn));
    nodeCount = ncL+ncR;
    if nn == 2
        nodeDiff = abs(ncL+ncR-nodeCount);
    end
    
    if options.debug
        hold off;
        plot(x(1:(N-idx+3)),uL/max(uL),'.-');
        hold on;
        plot(xflip(1:idx),uR/uR(end)*uL(end)/max(uL),'.-');
        pause(options.pauseTime);
    end
    match(nn) = duL(end)./uL(end-1)-duR(end)./uR(end-1);
    
    %%
    if abs(match(nn))<options.tol
        %If within tolerance, 
        solCount = solCount+1;
        Eout(solCount) = E(nn);
        u(:,solCount) = [uL(1:end-2);flip(uR(1:end-1)/uR(end-1)*uL(end-1))];
        break;
    end

    if nn ~= 1
        E(nn+1) = E(nn)-match(nn)*(E(nn)-E(nn-1))./(match(nn)-match(nn-1));
        err = abs(match(nn));
    end
    
    fprintf(1,'Iter: %02d, Error = %.3e, Energy = %.3f, Node Count = %d\n',nn,err,E(nn),nodeCount);
    nn = nn+1;
end
E = E(1:nn);
Eout = E(end);
u = [uL(1:end-2);flip(uR(1:end-1)/uR(end-1)*uL(end-1))];


end






