function u = makeSquareWellSolution(r,V0,r0,E)

if E>0
    k1 = sqrt(E-V0);
    u = zeros(size(r));
    u(r<=r0) = sin(k1*r(r<=r0));
    m = k1*cot(k1*r0);
    k = sqrt(E);
    K = (k*cos(k*r0)-m*sin(k*r0))./(m*cos(k*r0)+k*sin(k*r0));
    A = sin(k1*r0)./(sin(k*r0)+K*cos(k*r0));
    u(r>r0) = A*(sin(k*r(r>r0))+K*cos(k*r(r>r0)));
end




end