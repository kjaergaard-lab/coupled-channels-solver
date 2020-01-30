function idx=FindState(BV,vec)
% FindState Returns logical vector of size Nx1 where BV=vec where BV is NxM and vec is
% 1xM

idx=false(size(BV,1),1);
for nn=1:size(BV,1),
    idx(nn)=all(BV(nn,:)==vec);
end;

end