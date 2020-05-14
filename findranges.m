function Eout = findranges(x,Vfunc,Ein,recurse,ops,opt)

if nargin<6
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

opt.direction = 1;
opt.output = false;
opt.debug = false;

E = sort(Ein(:),'descend');
[~,~,nodes(1)] = manolopoulos(x,Vfunc,E(1),ops,opt);
[~,~,nodes(2)] = manolopoulos(x,Vfunc,E(2),ops,opt);
numBound = abs(nodes(1)-nodes(2));
if opt.debug
    fprintf(1,'E = [%.5f,%.5f], Bound states: %d\n',E(1),E(2),numBound);
end

Eout = [];
if numBound == 0
    Eout = [];
elseif numBound == 1 && recurse
    Enew = E;
    for nn=1:opt.rangeiter
        if opt.debug
            fprintf(1,'Iterating E = [%.5f,%.5f]\n',Enew(1),Enew(2));
        end
        Emid = (Enew(1)+Enew(2))/2;
        Etmp = findranges(x,Vfunc,[Enew(1),Emid],false,ops,opt);
        if isempty(Etmp)
            Etmp = findranges(x,Vfunc,[Emid,Enew(2)],false,ops,opt);
        end
        Enew = Etmp;
    end
    Eout = E;
elseif numBound == 1 && ~recurse
    Eout = E;
else
    Emid = (E(1)+E(2))/2;
    Enew = findranges(x,Vfunc,[E(1),Emid],recurse,ops,opt);
    if ~isempty(Enew)
        Eout = [Eout,Enew];
    end
    Enew = findranges(x,Vfunc,[Emid,E(2)],recurse,ops,opt);
    if ~isempty(Enew)
        Eout = [Eout,Enew];
    end
end


end