function V = RbPot(x,ell)

scale = const.mRb/const.hbar^2*1e-20*const.kb;
x = reshape(x,1,1,numel(x));
S = const.cm2K*RbBOPotentialsv2(x,0)*scale+ell*(ell+1)./(x).^2;
T = const.cm2K*RbBOPotentialsv2(x,1)*scale+ell*(ell+1)./(x).^2;

V = zeros(2,2,numel(x));
V(1,1,:) = S;
V(2,2,:) = T+15;
V(1,2,:) = 10e-2*(T-S);
% V(1,2,:) = 100*exp(-(x-5).^2./(2*2^2));
V(2,1,:) = V(1,2,:);


