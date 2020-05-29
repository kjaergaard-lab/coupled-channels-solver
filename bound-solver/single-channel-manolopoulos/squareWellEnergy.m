function squareWellEnergy(V0,x0)

k0 = sqrt(V0);
kE = linspace(0,k0,1e4);
k1 = sqrt(k0.^2-kE.^2);


figure(129);clf;
plot(kE.^2,k1+kE.*tan(k1.*x0),'.-');
ylim([-10,10]);


end