function y=RiccatiBessel(x,k,LMat,K,dopt)
% RiccatiBessel Calculates Riccati-Bessel functions and their derivatives
%   y=RiccatiBessel(x,k,LMat,K,dopt) x is position, k is the wavevector,
%   LMat is the matrix of L values, K=1,2 defines whether to use spherical
%   Hankel functions of type 1 (outgoing) or type 2 (incoming), and
%   dopt=0,1 calculates either the function (0) or its derivative (1) using
%   a combination of Bessel functions
%(dopt=1)
L=diag(LMat);
x=diag(k).*x;
if dopt==0,
    y=1i*sqrt(pi*x/2).*besselh(L+0.5,K,x);
else
    y=1i/2*sqrt(pi./(2.*x)).*(besselh(L+0.5,K,x)+x.*(besselh(L-0.5,K,x)-besselh(L+1.5,K,x)));
end;

if K==2,
    y=-y;
end;

y=diag(y);
y=sqrt(k)\y;

end