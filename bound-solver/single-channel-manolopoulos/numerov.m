function [u,du,x,idx,nodeCount] = numerov(x,ui,W,direction,options)

if nargin<5 || ~isfield(options,'stopAtRoot')
    options.stopAtRoot = 0;
end
if nargin<5 || ~isfield(options,'stopAtMax')
    options.stopAtMax = 0;
end
if nargin<5 || ~isfield(options,'maxTol')
    options.maxTol = 1e-3;
end
if nargin<5 || ~isfield(options,'stopXBelow')
    options.stopXBelow = max(x);
end

x = x(:);
W = W(:);
u = zeros(size(x));
if nargin>=4 && direction<0
    x = flip(x);
    W = flip(W);
end
dx = x(2)-x(1);
% dx = x(2:end)-x(1:end-1);

u(1:2) = ui(1:2);
nodeCount = 0;
for nn=3:numel(u)
    u(nn) = ((2+5/6*dx^2*W(nn-1))*u(nn-1)-(1-dx^2/12*W(nn-2))*u(nn-2))./(1-dx^2/12*W(nn));
%     h = x(nn-1)-x(nn-2);
%     c = (x(nn)-x(nn-1))/h;
%     h = dx(nn);
%     c = h./dx(nn-1);
%     u(nn) = ((1+c+h^2/12*(1+4*c+4*c^2+c^3)*W(nn-1))*u(nn-1)-(c-h^2/12*(c+c^2-c^3)*W(nn-2))*u(nn-2))/(1-h^2/12*(c+c^2-1)*W(nn));
    if sign(u(nn)) ~= sign(u(nn-1))
        nodeCount = nodeCount+1;
        if options.stopAtRoot
            break;
        end
    elseif options.stopAtMax && x(nn)<options.stopXBelow && isapprox((u(nn)-u(nn-1))/(u(nn)+u(nn-1)),0,options.maxTol)
        break;
    end
end
u = u(1:nn);
% dx = dx(1:(nn-1));
du = ((0.5+dx^2/12*W(3:nn)).*u(3:end)-(0.5+dx^2/12*W(1:nn-2)).*u(1:end-2))/dx;

% dy = u(2:end)-u(1:end-1);
% du = (dy(1:end-1).*dx(2:end).^2+dy(2:end).*dx(1:end-1).^2)./(dx(1:end-1).*dx(2:end).*(dx(1:end-1)+dx(2:end)));

% du = ((u(2:end-1)-u(1:end-2)).*(x(3:end)-x(2:end-1)).^2+(u(3:end)-u(2:end-1)).*(x(2:end-1)-x(1:end-2)).^2)...
%     ./((x(2:end-1)-x(1:end-2)).*(x(
x = x(1:nn);
idx = nn;








