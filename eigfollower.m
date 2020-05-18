function eigANew = eigfollower(r,Y)

Nch = size(Y,1);
[vecNew,eigNew] = eig(Y(:,:,1));
vecNew = gramschmidt(vecNew);
eigNew = diag(eigNew);
eigANew = zeros(Nch,size(Y,3));
eigANew(:,1) = eigNew;
vecANew = vecNew;

for nn=1:(size(Y,3)-1)
    vecAOld = vecANew;
    eigOld = eigNew;
    [vecNew,eigNew] = eig(Y(:,:,nn+1));
    vecNew = gramschmidt(vecNew);
    eigNew = diag(eigNew);
    for knew=1:Nch
        l2dist = 0;
        l2idx = 0;
        for kold=1:Nch
            l2Test = abs(vecNew(:,knew)'*vecAOld(:,kold));
            if l2Test>=l2dist
                l2dist = l2Test;
                l2idx = kold;
            end
        end
        vecANew(:,l2idx) = vecNew(:,knew);
        eigANew(l2idx,nn+1) = eigNew(knew);
    end
end