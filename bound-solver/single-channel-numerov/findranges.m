function Eout = findranges(x,V,Ein,options)

if nargin<4 || ~isfield(options,'recurse')
    options.recurse = true;
end
if nargin<4 || ~isfield(options,'iter')
    options.iter = 4;
end

Ein = Ein(:);
[~,~,~,~,nodes(1)] = numerov(x,[0,1e-20],V-Ein(1),1);
[~,~,~,~,nodes(2)] = numerov(x,[0,1e-20],V-Ein(2),1);
numBound = abs(nodes(1)-nodes(2));

Eout = [];
if numBound == 0
    Eout = [];
elseif numBound == 1 && options.recurse
    options2 = options;
    options2.recurse = false;
    for nn=1:options.iter
        Emid = (Ein(1)+Ein(2))/2;
        Enew = findranges(x,V,[Ein(1),Emid],options2);
        if ~isempty(Enew)
            Ein = Enew;
        else
            Ein = findranges(x,V,[Emid,Ein(2)],options2);
        end
    end
    Eout = Ein;
elseif numBound == 1 && ~options.recurse
    Eout = Ein;
else
    Emid = (Ein(1)+Ein(2))/2;
    Enew = findranges(x,V,[Ein(1),Emid],options);
    if ~isempty(Enew)
        Eout = [Eout,Enew];
    end
    Enew = findranges(x,V,[Emid,Ein(2)],options);
    if ~isempty(Enew)
        Eout = [Eout,Enew];
    end
end


end