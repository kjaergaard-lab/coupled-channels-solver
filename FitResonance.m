function [params,fitfun]=FitResonance(x,y,sel)
% FitResonance Fits a resonance to four different options
%   [params,fitfun]=FitResonance(x,y,sel) with x,y the x and y values and
%   sel=0,1,2,3,4 selecting the fit type.  Sel=0: fits phase to arctangent.
%   Sel=1 Fits cross-section (or abs(T).^2) to Beutler-Fano profile. Sel=2
%   fits to Lorentzian. Sel=3 fits real part of scattering length to appropriate form
%   for magnetic Feshbach resonance with inelastic decay. Sel=4 fits
%   scattering length to form for no inelastic decay.
x=x(:);
y=y(:);
% xplot=sort(unique([x',linspace(min(x),max(x),1e4)]))';
xplot=linspace(min(x),max(x),1e4);
if sel==0,
    func=@(c,x) unwrap((c(1)+atan(0.5*c(2)./(x-c(3))))*2)/2;
    lb=[-pi,-100,min(x)];ub=[pi,100,max(x)];
    [~,idx]=max(abs(diff(y)));
    guess=[y(1),0.05,x(idx)];
    options=optimset('Display','off');
    params=lsqcurvefit(func,guess,x,y,lb,ub,options);
    fitfun=func(params,xplot);
elseif sel==1,
    func=@(c,x) c(1)*sin(c(2)+atan(0.5*c(3)./(x-c(4)))).^2;
    
    [m,idx]=max(y);
%     m=4;
    [~,idx2]=min(y);
    d=sign(x(idx)-x(idx2))*asin(sqrt((y(1)/4)));
    B0=cos(d).^2*x(idx)+sin(d).^2*x(idx2);
    G=sin(2*d)*(x(idx)-x(idx2));
    x = x-B0;
    guess=[m,d,G,B0*0];
    lb=[0,-2*pi,0,min(x)];ub=[4,2*pi,150,max(x)];
    options=optimset('Display','off','TolX',1e-11,'TolFun',1e-11,'MaxFunEvals',1e5,'Maxiter',1e4);
    
    params=lsqcurvefit(func,guess,x,y,lb,ub,options);
    params(4) = params(4)+B0;
    x = x+B0;
    fitfun=func(params,xplot);
elseif sel==2,
    func=@(c,x) c(1)*c(2).^2./4./(c(2).^2./4+(x-c(3)).^2);
    lb=[0,0,min(x)];ub=[2e3*max(y),1,max(x)];
    [m,idx]=max(abs(y));
    guess=[m,0.1,x(idx)];
    options=optimset('Display','off');
    params=lsqcurvefit(func,guess,x,y,lb,ub,options);
    fitfun=func(params,xplot);
elseif sel==3,
    func=@(c,x) c(1).*(1-c(2).*(x-c(3))./((x-c(3)).^2+c(4).^2/4));
    lb=[0,0,min(x),0];ub=[1e5,2000,max(x),1];
%     lb=[0.5,-pi,0,min(x)];ub=[1,pi,1,max(x)];
%     [m,idx]=max(abs(y));
    [~,idx]=max(diff(y));
    guess=[y(1),0.005,x(idx),0.0005];
    options=optimset('Display','off','TolX',1e-9,'TolFun',1e-9);
    params=lsqcurvefit(func,guess,x,y,lb,ub,options);
    fitfun=func(params,xplot);
elseif sel==4,
    func=@(c,x) c(1).*(1-c(2)./(x-c(3)));
    lb=[0,0,min(x)];ub=[1e5,20,max(x)];
%     lb=[0.5,-pi,0,min(x)];ub=[1,pi,1,max(x)];
    [m,idx]=max(abs(y));
%     [~,idx]=max(diff(y));
    guess=[y(1),0.005,x(idx)];
    options=optimset('Display','off','TolX',1e-9,'TolFun',1e-9);
    params=lsqcurvefit(func,guess,x,y,lb,ub,options);
    fitfun=func(params,xplot);
end;

figure(20);clf(20);
R=func(params,x)-y;
% subplot(3,1,1:2);
plot(x,y,'b.');
hold on;
plot(xplot,fitfun,'r--');
% title(num2str(sum(R.^2)));
hold off;
% subplot(3,1,3);
% plot(x,R,'b.-');

