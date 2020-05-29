function [S,BVSymInt,ia] = SymmetrizeSMatrixArbL(Sin,bvint,symmetry)
% SymmetrizeSMatrixArbL symmetrizes S-matrix for identical particles
% [S,BVSymInt]=SymmetrizeSMatrixArbL(Sin,bvint,symmetry) Sin is input S
% matrix, BVint is input interal basis vectors, and symmetry=1,-1 defines
% whether it is symmetric or anti-symmetric symmetrization.  S is output S
% matrix and BVSymInt is the compressed form of BVint

Nchannels = size(Sin,1);
[~,idx] = sort(bvint(:,5:6),2);
tmp = zeros(size(bvint));
for nn=1:Nchannels
    tt1 = [3,4];tt2 = [5,6];
    tmp(nn,:) = [bvint(nn,1:2) bvint(nn,tt1(idx(nn,:))) bvint(nn,tt2(idx(nn,:)))];
end
[BVSymInt,ia] = unique(tmp,'stable','rows');
NSymStates = size(BVSymInt,1);
S = zeros(NSymStates,NSymStates,size(Sin,3));

for row=1:NSymStates
    for col=1:NSymStates
        idxR1 = all(repmat(BVSymInt(row,:),Nchannels,1)==bvint,2);
        idxR2 = all(repmat(BVSymInt(row,[1:2,4,3,6,5]),Nchannels,1)==bvint,2);
        idxC1 = all(repmat(BVSymInt(col,:),Nchannels,1)==bvint,2);
        idxC2 = all(repmat(BVSymInt(col,[1:2,4,3,6,5]),Nchannels,1)==bvint,2);
        S(row,col,:) = Sin(idxR1,idxC1,:)+symmetry*(-1).^(BVSymInt(col,1))*Sin(idxR1,idxC2,:)+symmetry*(-1).^(BVSymInt(row,1))*Sin(idxR2,idxC1,:)+(-1).^(BVSymInt(row,1)+BVSymInt(col,1))*Sin(idxR2,idxC2,:);
        dR = BVSymInt(row,5)==BVSymInt(row,6);
        dC = BVSymInt(col,5)==BVSymInt(col,6);
        S(row,col,:) = S(row,col,:)./(sqrt(4*(1+dR).*(1+dC)));
        S(row,col,:) = S(row,col,:)+(S(row,col,:)==0 & row==col);
    end
end
