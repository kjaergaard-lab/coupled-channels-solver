function Eout = findranges(x,Vfunc,Ein,recurse,iter,opt)

if nargin<6
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

opt.direction = 1;

Ein = sort(Ein(:),'descend');
[~,~,nodes(1)] = manolopoulos(x,Vfunc,Ein(1),opt);
[~,~,nodes(2)] = manolopoulos(x,Vfunc,Ein(2),opt);
numBound = abs(nodes(1)-nodes(2));

Eout = [];
if numBound == 0
    Eout = [];
elseif numBound == 1 && recurse
    for nn=1:iter
        Emid = (Ein(1)+Ein(2))/2;
        Enew = findranges(x,Vfunc,[Ein(1),Emid],false,iter,opt);
        if ~isempty(Enew)
            Ein = Enew;
        else
            Ein = findranges(x,Vfunc,[Emid,Ein(2)],false,iter,opt);
        end
    end
    Eout = Ein;
elseif numBound == 1 && ~recurse
    Eout = Ein;
else
    Emid = (Ein(1)+Ein(2))/2;
    Enew = findranges(x,Vfunc,[Ein(1),Emid],recurse,iter,opt);
    if ~isempty(Enew)
        Eout = [Eout,Enew];
    end
    Enew = findranges(x,Vfunc,[Emid,Ein(2)],recurse,iter,opt);
    if ~isempty(Enew)
        Eout = [Eout,Enew];
    end
end


end