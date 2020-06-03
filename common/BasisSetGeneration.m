function basis = BasisSetGeneration(AtomPair,LVec,filename)
% BasisSetGeneration Creates BasisSetFile for use in MultiChannel
%   BasisSetGeneration(AtomPair,LVec,filename)
%   AtomPair is a string and one of 'RbRb', 'KK', or 'KRb'.
%   LVec defines the range of angular momentum quantum numbers to include,
%   as for instance [0:2], or [1,3,5]
%   filename is the name of the file that is saved.  Default is [AtomPair,
%   ' Basis Set L=LVec']
%   See also MultiChannel, ClebschGordan

%% Preamble
if nargin==2
    Lstr = ProcessLstr(LVec);
    filename = [AtomPair,' Basis Set ',Lstr];
end

basis = atompairbasis(AtomPair,LVec);
save(filename,'basis');

end

function Lstr = ProcessLstr(LVec)
% Produces a string from LVec that compresses concurrent numbers
% So LVec=[0:2,5,8:10] becomes 'L=0-2,5,8-10'

LVec = sort(unique(LVec));
if numel(LVec)==2
    if diff(LVec)==1
        Lstr = ['L=',num2str(LVec(1)),'-',num2str(LVec(2))];
    else
        Lstr = ['L=',num2str(LVec(1)),',',num2str(LVec(2))];
    end
else
    Lstr = ['L=',num2str(LVec(1))];
end

for nn=2:numel(LVec)-1
    dL = diff(LVec(nn-1:nn+1));
    if dL(1)==1 && dL(2)~=1
        Lstr = [Lstr,'-',num2str(LVec(nn))]; %#ok<*AGROW>
    elseif dL(1)~=1 && dL(2)==1
        Lstr = [Lstr,',',num2str(LVec(nn))];
    elseif dL(1)~=1 && dL(2)~=1
        Lstr = [Lstr,',',num2str(LVec(nn))];
    end
    
    if nn==numel(LVec)-1
        if dL(2)==1
           Lstr = [Lstr,'-',num2str(LVec(nn+1))];
        else
            Lstr = [Lstr,',',num2str(LVec(nn+1))];
        end           
    end
end
end