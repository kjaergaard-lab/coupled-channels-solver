function varargout = manolopoulos(r,Vfunc,E,varargin)

%% Default values
opt.stopAtRoot = false;
opt.stopAtR = false;
opt.stopAfterR = false;
opt.stopR = 5;
opt.direction = 1;
opt.makeOutput = false;

%% Parse arguments
if numel(varargin)==1 && isa(varargin{1},'struct')
    options = flattenstruct(varargin{1});
elseif numel(varargin)>0 && mod(numel(varargin),2)~=0
    error('Must supply name/value pairs for additional arguments');
else
    options = varargin;
end

for nn=1:2:numel(options)
    if ~ischar(options{nn})
        error('Optional argument %d is not a string',nn);
    else
        v = options{nn+1};
        switch lower(options{nn})
            case 'stopatroot'
                opt.stopAtRoot = v;
            case 'stopatr'
                opt.stopAtR = v;
            case 'stopafterr'
                opt.stopAfterR = v;
            case 'stopr'
                opt.stopR = v;
            case 'direction'
                opt.direction = v;
            case 'makeoutput'
                opt.makeOutput = v;
%             otherwise
%                 error('Option %s not known',options{nn});
        end
    end
end

if opt.direction>0
    Ynew = 1e20;
    unew = 1e-20;
elseif opt.direction<0
    r = flip(r);
    Ynew = -sqrt(Vfunc(r(1))-E);
    unew = 1e-20;
end
Yold = Ynew;

if opt.makeOutput
    Y = zeros(numel(r),1);
    Y(1) = Ynew;
    u = zeros(numel(r),1);
    u(1) = unew;
end

nodes = 0;
dr = diff(r);
h = dr/2;
M = Vfunc(r)-E;
M2 = Vfunc(r(1:end-1)+h)-E;
for nn=1:numel(r)-1
%     dr = r(nn+1)-r(nn);
%     h = dr/2;
    p2 = M2(nn);
    p = sqrt(abs(p2));
    if p2 == 0
        y1 = 1/h(nn);
        y2 = 1/h(nn);
    else
        y1 = p.*coth(p*h(nn)).*(p2>0)+p.*cot(p*h(nn)).*(p2<0);
        y2 = p.*csch(p*h(nn)).*(p2>0)+p.*csc(p*h(nn)).*(p2<0);
    end
%     y1 = 1/h*(cothc(p*h).*(p2>0)+cotc(p*h).*(p2<0));
%     y2 = 1/h*(cschc(p*h).*(p2>0)+cscc(p*h).*(p2<0));
    
    Qa = h(nn)/3*(M(nn)-p2);
%         Qc = 4/(h*(1-h^2/6*(M2(nn)-p2)))-4/h;
    Qc = 0;
    Qb = h(nn)/3*(M(nn+1)-p2);

    Yold2 = Yold;
    Yold = Ynew;
    Ytmp = (y1+Qc)-y2/(Ynew+y1+Qa)*y2;
    utmp = (Ynew+y1+Qa)/y2*unew;
    Ynew = (y1+Qb)-y2/(Ytmp+y1+Qc)*y2;
    unew = (Ytmp+y1+Qc)/y2*utmp;
    if makeOutput
        Y(nn+1) = Ynew;
        u(nn+1) = unew;
    end
    
    if sign(Yold) ~= sign(Ynew)
        if sign(Ynew-Yold) ~= sign(Yold-Yold2)
            nodes = nodes+1;
        elseif opt.stopAtRoot
            if opt.direction>0 && r(nn+1)>opt.stopR
                break;
            elseif opt.direction<0 && r(nn+1)<opt.stopR
                break;
            end
        end
    end

    if opt.stopAtR && isapprox(r(nn+1),opt.stopR,1e-10)
        break;
    elseif opt.stopAfterR
        if opt.direction>0 && r(nn+1)>=opt.stopR
            break;
        elseif opt.direction<0 && r(nn+1)<=opt.stopR
            break;
        end
    end
end

if opt.makeOutput
    Y = Y(1:(nn+1));
    u = u(1:(nn+1));
end


varargout{1} = Ynew;
varargout{2} = r(nn+1);
varargout{3} = nodes;
if opt.makeOutput && nargout>3
    varargout{4} = r(1:(nn+1));
    varargout{5} = Y;
    varargout{6} = u;
end

end
    
