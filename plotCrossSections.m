clear cs str

states = unique(sm.BV(:,3:4),'rows');

jj = 1;
clf;
for nn=1:size(states,1)
    tmp = CalculatePartialCrossSection(Tsym,sm.BV,sm.k,InitStateLabel(3:4),states(nn,1:2));
    if tmp(end)~=0
        cs(:,jj) = tmp*(1+(states(nn,1)==states(nn,2)));
        str{jj} = sprintf('%d-%d',states(nn,1),states(nn,2));
        plot(sm.E,cs(:,jj),'.-');
        hold on;
        jj = jj+1;
    end
end
hold off;
legend(str);




