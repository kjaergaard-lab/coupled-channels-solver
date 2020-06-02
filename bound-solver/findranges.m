function Eout = findranges(r,Ein,recurse,ops,opt)
% FINDRANGES Determines energy ranges containing exactly one bound state
%
%   Eout = FINDRANGES(r,Ein,recurse,ops,opt) uses the vector of positions
%   r to find all energy ranges Eout between which there is exactly one
%   bound state based on node counting.  Eout is a 2xNbound vector.
%
%   Set recurse = true in the initial call to this function so that it
%   searches for all bound state ranges.
%
%   OPS is an instance of ALKALIOPERATORS
%
%   OPT is an instance of BOUNDOPTIONS.  Relevant boundoptions properties
%   are:
%
%   rangeiter: number of iterations to use for narrowing energy ranges.
%   More iterations means the bound state solver is typically more likely
%   to converge on a solution without using fsolve(), but it takes longer

if nargin<5
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

opt.direction = 1;
opt.output = false;
opt.debug = false;

E = sort(Ein(:),'descend');
[~,~,nodes(1)] = manolopoulos_bound(r,E(1),ops,opt);
[~,~,nodes(2)] = manolopoulos_bound(r,E(2),ops,opt);
numBound = abs(nodes(1)-nodes(2));
if opt.debug
    fprintf(1,'E = [%.5f,%.5f], Bound states: %d\n',E(1),E(2),numBound);
end

Eout = [];
if numBound == 0
    Eout = [];
elseif numBound == 1 && recurse
    for nn=1:opt.rangeiter
        if opt.debug
            fprintf(1,'Iterating E = [%.5f,%.5f]\n',E(1),E(2));
        end
        Emid = (E(1)+E(2))/2;
        Etmp = findranges(r,[E(1),Emid],false,ops,opt);
        if isempty(Etmp)
            Etmp = findranges(r,[Emid,E(2)],false,ops,opt);
        end
        E = Etmp;
    end
    Eout = E;
elseif numBound == 1 && ~recurse
    Eout = E;
else
    Emid = (E(1)+E(2))/2;
    Enew = findranges(r,[E(1),Emid],recurse,ops,opt);
    if ~isempty(Enew)
        Eout = [Eout,Enew];
    end
    Enew = findranges(r,[Emid,E(2)],recurse,ops,opt);
    if ~isempty(Enew)
        Eout = [Eout,Enew];
    end
end


end