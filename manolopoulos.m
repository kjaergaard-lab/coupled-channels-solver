function varargout = manolopoulos(r,Vfunc,E,opt)

%% Parse arguments
if nargin<4
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

%% Prep
I = eye(size(Vfunc(r(1))));
Nch = size(I,1);
if opt.direction>0
    Ynew = sqrt(Vfunc(r(1))-E*I);
elseif opt.direction<0
    r = flip(r);
    Ynew = -sqrt(Vfunc(r(1))-E);
end
% Yold = Ynew;
eigNew = eig(Ynew);
negEigNew = sum(eigNew<0);
detNew = prod(eigNew);
eigOld = eigNew;
negEigOld = negEigNew;
detOld = detNew;
[~,minIndex] = min(abs(eigNew));
minEigNew = eigNew(minIndex);
[~,maxIndex] = max(abs(eigNew));
maxEigNew = eigNew(maxIndex);
maxEigOld = maxEigNew;

% detOld = det(Ynew);
% detNew = detOld;

if opt.output
    Y = zeros(Nch,Nch,numel(r));
    Y(:,:,1) = Ynew;
    Zout = zeros(Nch,Nch,numel(r));
end

nodes = zeros(numel(r),1);
antinodes = zeros(numel(r),1);
dbg.test = zeros(numel(r),1);
dr = diff(r);
h = dr/2;
Epot = repmat(E*I,1,1,numel(r));
M = Vfunc(r)-Epot;
M2 = Vfunc(r(1:end-1)+h)-Epot(:,:,1:end-1);

%% Solve
for nn=1:numel(r)-1
    p2 = diag(M2(:,:,nn));
    p = sqrt(abs(p2));
    y1 = p.*coth(p*h(nn)).*(p2>0)+p.*cot(p*h(nn)).*(p2<0)+1./h(nn).*(p2==0);
    y2 = p.*csch(p*h(nn)).*(p2>0)+p.*csc(p*h(nn)).*(p2<0)+1./h(nn).*(p2==0);
    
    y1 = diag(y1);
    y2 = diag(y2);
    Mref = diag(p2);
    
    Qa = h(nn)/3*(M(:,:,nn)-Mref);
    Qc = 4/h(nn)*((I-h(nn)^2/6*(M2(:,:,nn)-Mref))\I)-4/h(nn)*I;
    Qb = h(nn)/3*(M(:,:,nn+1)-Mref);

    eigOld2 = eigOld;
    negEigOld2 = negEigOld;
    detOld2 = detOld;
    eigOld = eigNew;
    negEigOld = negEigNew;
    detOld = detNew;
    minEigOld = minEigNew;
    maxEigOld2 = maxEigOld;
    maxEigOld = maxEigNew;

    Zac = (Ynew+y1+Qa)\y2;
    Yc = (y1+Qc)-y2*Zac;
    Zcb = (Yc+y1+Qc)\y2;
    Ynew = (y1+Qb)-y2*Zcb;

    eigNew = eig(Ynew);
    negEigNew = sum(eigNew<0);
    detNew = prod(eigNew);
%     [~,minIndex] = min(abs(eigNew));
    minEig = Inf;
    maxEig = 0;
    minIndex = 0;
%     maxIndex = 0;
    for kk=1:numel(eigNew)
        eigTest = abs(eigNew(kk));
        if eigTest<minEig
            minEig = eigTest;
            minIndex = kk;
        elseif eigTest>maxEig
            maxEig = eigTest;
%             maxIndex = kk;
        end
    end
        
    minEigNew = eigNew(minIndex);
%     [~,maxIndex] = max(abs(eigNew));
%     maxEigNew = eigNew(maxIndex);
    
    if opt.output
        Y(:,:,nn+1) = Ynew;
        Zout(:,:,nn+1) = Zac*Zcb;
    end
    
    nodes(nn+1) = nodes(nn);
    antinodes(nn+1) = antinodes(nn);
    if (negEigNew>negEigOld && opt.direction>0) || (negEigNew<negEigOld && opt.direction<0)
        if opt.stopAtRoot
            if opt.direction>0 && r(nn+1)>opt.stopR
                break;
            elseif opt.direction<0 && r(nn+1)<opt.stopR
                break;
            end
        end
    elseif (negEigNew<negEigOld && opt.direction>0) || (negEigNew>negEigOld && opt.direction<0)
        nodes(nn+1) = nodes(nn)+opt.direction*(negEigOld-negEigNew);
    elseif (minEigNew<0 && minEigOld>0 && abs(maxEig/minEigNew)>10 && opt.direction>0) || (minEigNew>0 && minEigOld<0 && opt.direction<0)
        nodes(nn+1) = nodes(nn)+1;
    end
    
    if (negEigNew>negEigOld && opt.direction>0) || (negEigNew<negEigOld && opt.direction<0)
        antinodes(nn+1) = antinodes(nn)+opt.direction*(negEigNew-negEigOld);
    elseif (negEigNew==negEigOld && minEigNew<0 && minEigOld>0 && abs(maxEig/minEigNew)>10 && opt.direction>0) || (negEigNew==negEigOld && minEigNew>0 && minEigOld<0 && opt.direction<0)
        antinodes(nn+1) = antinodes(nn)+1;
    end
    

    if opt.stopAtR && abs(r(nn+1)-opt.stopR)<1e-10
        break;
    elseif opt.stopAfterR
        if opt.direction>0 && r(nn+1)>=opt.stopR
            break;
        elseif opt.direction<0 && r(nn+1)<=opt.stopR
            break;
        end
    end
end

if opt.output
    Y = Y(:,:,1:(nn+1));
    Zout = Zout(:,:,1:(nn+1));
end

dbg.nodes = nodes;
dbg.antinodes = antinodes;
nodes = nodes(nn+1);

varargout{1} = Ynew;
varargout{2} = r(nn+1);
varargout{3} = nodes;
if opt.output && nargout>3
    varargout{4} = r(1:(nn+1));
    varargout{5} = Y;
    varargout{6} = Zout;
    varargout{7} = dbg;
end

end
    
