function eigANew = eigfollower(r,Y)

Nch = size(Y,1);
[vecNew,eigNew] = eig(Y(:,:,1));
eigNew = diag(eigNew);
eigANew = zeros(Nch,size(Y,3));
eigANew(:,1) = eigNew;
vecANew = vecNew;

for nn=1:(size(Y,3)-1)
    vecAOld = vecANew;
    eigOld = eigNew;
    [vecNew,eigNew] = eig(Y(:,:,nn+1));
    eigNew = diag(eigNew);
%     eigAOld = eigANew;
%     eigANew = eigNew;
    for knew=1:Nch
        l2dist = 0;
        l2idx = 0;
        for kold=1:Nch
%             l2Test = 0;
%             for tt=1:Nch
%                 l2Test = l2Test + (vecNew(tt,knew)-vecAOld(tt,kold))^2;
%             end
%             l2Test = norm(vecNew(:,knew)-vecAOld(:,kold),2);
            l2Test = abs(vecNew(:,knew)'*vecAOld(:,kold));
            if l2Test>l2dist
                l2dist = l2Test;
                l2idx = kold;
            end
        end
        vecANew(:,knew) = vecNew(:,l2idx);
        eigANew(knew,nn+1) = eigNew(l2idx);
    end
end