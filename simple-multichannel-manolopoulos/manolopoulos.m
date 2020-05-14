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
    opt.blocks = rot90(opt.blocks(end,end)+1-opt.blocks,2);
    Ynew = -sqrt(Vfunc(r(1))-E*I);
    
end
Ynew = diag(diag(Ynew));

[vecNew,eigNew] = eig(Ynew);
eigNew = diag(eigNew);
negEigNew = sum(eigNew<0);
eigAOld = eigNew;
eigANew = eigNew;
vecANew = vecNew;


if opt.output
    Y = zeros(Nch,Nch,numel(r));
    Y(:,:,1) = Ynew;
    Zout = zeros(Nch,Nch,numel(r));
end

breakFlag = false;
nodes = zeros(numel(r),1);
dbg.test = zeros(numel(r),1);
gcnt = 1;

for bb=1:size(opt.blocks,1)
    Epot = repmat(E*I,1,1,opt.blocks(bb,2)-opt.blocks(bb,1)+1);
    rb = r(opt.blocks(bb,1):opt.blocks(bb,2));
    h = (rb(2)-rb(1))/2;
    M = Vfunc(rb)-Epot;
%     M2 = Vfunc(rb(1:end-1)+h)-Epot(:,:,1:end-1);
    M2 = Vfunc(rb+h)-Epot;

    %% Solve
    for nn=1:numel(rb)-1
        p2 = diag(M2(:,:,nn));
        p = sqrt(abs(p2));
        y1 = p.*coth(p*h).*(p2>0)+p.*cot(p*h).*(p2<0)+1./h.*(p2==0);
        y2 = p.*csch(p*h).*(p2>0)+p.*csc(p*h).*(p2<0)+1./h.*(p2==0);
        
        y1 = diag(y1);
        y2 = diag(y2);
        Mref = diag(p2);
        
        Qa = h/3*(M(:,:,nn)-Mref);
        Qc = 4/h*((I-h^2/6*(M2(:,:,nn)-Mref))\I)-4/h*I;
        Qb = h/3*(M(:,:,nn+1)-Mref);
        
        Zac = (Ynew+y1+Qa)\y2;
        Yc = (y1+Qc)-y2*Zac;
        Zcb = (Yc+y1+Qc)\y2;
        Ynew = (y1+Qb)-y2*Zcb;
        
        vecAOld = vecANew;
        eigAOld2 = eigAOld;
        eigAOld = eigANew;
        [vecNew,eigNew] = eig(Ynew);
        eigNew = diag(eigNew);
        minEigNew = Inf;
        maxEigNew = 0;
        for kk=1:Nch
            eigTest = abs(eigNew(kk));
            if eigTest<minEigNew
                minEigNew = eigTest;
            elseif eigTest>maxEigNew
                maxEigNew = eigTest;
            end
        end
        
        if opt.output
            Y(:,:,gcnt+1) = Ynew;
            Zout(:,:,gcnt+1) = Zac*Zcb;
        end
        
        for knew=1:Nch
            l2dist = 0;
            l2idx = 0;
            for kold=1:Nch
                l2Test = abs(vecNew(:,knew)'*vecAOld(:,kold));
                if l2Test>=l2dist
                    l2dist = l2Test;
                    l2idx = kold;
                end
            end
            vecANew(:,knew) = vecNew(:,l2idx);
            eigANew(knew) = eigNew(l2idx);
        end
        
        dnodes = 0;
        for kk=1:Nch
            if sign(eigAOld(kk)) ~= sign(eigANew(kk))
                if (sign(eigANew(kk)-eigAOld(kk)) ~= sign(eigAOld(kk)-eigAOld2(kk))) && (abs(eigANew(kk))>0.1 || abs(eigAOld(kk))>0.1)
                    dnodes = dnodes+1;
                elseif opt.stopAtRoot
                    if opt.direction>0 && rb(nn+1)>opt.stopR
                        breakFlag = true;
                    elseif opt.direction<0 && rb(nn+1)<opt.stopR
                        breakFlag = true;
                    end
                end
            end
        end
        nodes(gcnt+1) = nodes(gcnt)+dnodes;
        
        if opt.stopAtR && abs(rb(nn+1)-opt.stopR)<1e-10
            breakFlag = true;
        elseif opt.stopAfterR
            if opt.direction>0 && rb(nn+1)>=opt.stopR
                breakFlag = true;
            elseif opt.direction<0 && rb(nn+1)<=opt.stopR
                breakFlag = true;
            end
        end
        
        gcnt = gcnt+1;
        
        if breakFlag
            break;
        end
    end
    
    if breakFlag
        break;
    end
    
end

if opt.output
    Y = Y(:,:,1:gcnt);
    Zout = Zout(:,:,1:gcnt);
end

dbg.nodes = nodes;
nodes = nodes(gcnt);

varargout{1} = Ynew;
varargout{2} = r(gcnt);
varargout{3} = nodes;
if opt.output && nargout>3
    varargout{4} = r(1:gcnt);
    varargout{5} = Y;
    varargout{6} = Zout;
    varargout{7} = dbg;
end

end
    
