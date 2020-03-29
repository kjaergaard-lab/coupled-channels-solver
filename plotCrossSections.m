clear cs str intState

[~,ia] = unique(BVSymInt(:,3:4),'rows');
intState = BVSymInt(ia,:);
jj = 1;
clf;
for nn=1:size(intState,1)
    tmp = CalculatePartialCrossSection(Tsym,BVSymInt,sm.k,InitStateLabel(3:4),intState(nn,3:4));
    if tmp(end)~=0
        cs(:,jj) = tmp*(1+1*(intState(nn,3)==intState(nn,4)));
        plot(sm.E,cs(:,jj),'.-');
        hold on;
        str{jj} = sprintf('%d-%d',intState(nn,3:4));
        jj=jj+1;
    end
end
hold off;
legend(str);

