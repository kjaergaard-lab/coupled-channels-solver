function [S,BVSymInt]=SymmetrizeSMatrixArbL(Sin,BVint,SymType)
% SymmetrizeSMatrixArbL symmetrizes S-matrix for identical particles
% [S,BVSymInt]=SymmetrizeSMatrixArbL(Sin,BVint,SymType) Sin is input S
% matrix, BVint is input interal basis vectors, and SymType=1,-1 defines
% whether it is symmetric or anti-symmetric symmetrization.  S is output S
% matrix and BVSymInt is the compressed form of BVint

NChannels=size(Sin,1);
[~,idx]=sort(BVint(:,3:4),2);
tmp=zeros(size(BVint));
for nn=1:NChannels
    tmp(nn,:)=[BVint(nn,1:2) BVint(nn,(idx(nn,:)==1)*3+(idx(nn,:)==2)*4)];
end
BVSymInt=unique(tmp,'rows');
NSymStates=size(BVSymInt,1);
S=zeros(NSymStates,NSymStates,size(Sin,3));

for row=1:NSymStates
    for col=1:NSymStates
        idxR1=all(repmat(BVSymInt(row,:),NChannels,1)==BVint,2);
        idxR2=all(repmat(BVSymInt(row,[1:2,4,3]),NChannels,1)==BVint,2);
        idxC1=all(repmat(BVSymInt(col,:),NChannels,1)==BVint,2);
        idxC2=all(repmat(BVSymInt(col,[1:2,4,3]),NChannels,1)==BVint,2);
        S(row,col,:)=Sin(idxR1,idxC1,:)+SymType*(-1).^(BVSymInt(col,1))*Sin(idxR1,idxC2,:)+SymType*(-1).^(BVSymInt(row,1))*Sin(idxR2,idxC1,:)+(-1).^(BVSymInt(row,1)+BVSymInt(col,1))*Sin(idxR2,idxC2,:);
        dR=BVSymInt(row,3)==BVSymInt(row,4);
        dC=BVSymInt(col,3)==BVSymInt(col,4);
        S(row,col,:)=S(row,col,:)./(sqrt(4*(1+dR).*(1+dC)));
        S(row,col,:)=S(row,col,:)+(S(row,col,:)==0 & row==col);
    end
end
