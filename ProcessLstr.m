function Lstr=ProcessLstr(LVec)
% LVec=[0,2,4:6,7 10];

LVec=sort(unique(LVec));

if numel(LVec)==2,
    if diff(LVec)==1,
        Lstr=['L=',num2str(LVec(1)),'-',num2str(LVec(2))];
    else
        Lstr=['L=',num2str(LVec(1)),',',num2str(LVec(2))];
    end;
else
    Lstr=['L=',num2str(LVec(1))];
end;

for nn=2:numel(LVec)-1,
    dL=diff(LVec(nn-1:nn+1));
    if dL(1)==1 && dL(2)~=1
        Lstr=[Lstr,'-',num2str(LVec(nn))];
    elseif dL(1)~=1 && dL(2)==1,
        Lstr=[Lstr,',',num2str(LVec(nn))];
    elseif dL(1)~=1 && dL(2)~=1,
        Lstr=[Lstr,',',num2str(LVec(nn))];
    end;
    
    if nn==numel(LVec)-1,
        if dL(2)==1,
           Lstr=[Lstr,'-',num2str(LVec(nn+1))];
        else
            Lstr=[Lstr,',',num2str(LVec(nn+1))];
        end;            
    end;
end;


