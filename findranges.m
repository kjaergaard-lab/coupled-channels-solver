function Eout = findranges(x,Vfunc,Ein,options)

%% Default values
if nargin<4 || ~isfield(options,'recurse')
    options.recurse = true;
end
if nargin<4 || ~isfield(options,'iter')
    options.iter = 4;
end
options.direction = 1;

Ein = Ein(:);
[~,~,nodes(1)] = manolopoulos(x,Vfunc,Ein(1),options);
[~,~,nodes(2)] = manolopoulos(x,Vfunc,Ein(2),options);
numBound = abs(nodes(1)-nodes(2));

Eout = [];
if numBound == 0
    Eout = [];
elseif numBound == 1 && options.recurse
    options2 = options;
    options2.recurse = false;
    for nn=1:options.iter
        Emid = (Ein(1)+Ein(2))/2;
        Enew = findranges(x,Vfunc,[Ein(1),Emid],options2);
        if ~isempty(Enew)
            Ein = Enew;
        else
            Ein = findranges(x,Vfunc,[Emid,Ein(2)],options2);
        end
    end
    Eout = Ein;
elseif numBound == 1 && ~options.recurse
    Eout = Ein;
else
    Emid = (Ein(1)+Ein(2))/2;
    Enew = findranges(x,Vfunc,[Ein(1),Emid],options);
    if ~isempty(Enew)
        Eout = [Eout,Enew];
    end
    Enew = findranges(x,Vfunc,[Emid,Ein(2)],options);
    if ~isempty(Enew)
        Eout = [Eout,Enew];
    end
end


end