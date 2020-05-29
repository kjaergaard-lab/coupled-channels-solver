function params=LoopFit(B,y,sel)
%Fits data in a loop
N=size(y,2);

for nn=1:N,
%     params(nn,:)=FitResonance(B,y(:,nn)./max(y(:,nn)),1);
    params(nn,:)=FitResonance(B,y(:,nn),sel);
    pause(0.05);
end;



