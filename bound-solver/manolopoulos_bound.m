function varargout = manolopoulos_bound(r,E,ops,opt)
%manolopoulos_bound Integrates the matrix Riccati equation using the
%improved log-derivative method of D. E. Manolopoulos for bound state
%problems
%
%   [Y,rf,nodes] = manolopoulos_bound(r,E,ops,opt) integrates the matrix Riccati
%   equation outwards assuming an energy E.
%
%   ops is an instance of ALKALIOPERATORS
%
%   opt is an instance of SCATTOPTIONS.  Relevant properties are:
%
%   output: set to true to generate wavefunctions
%   direction: direction to integrate.  Set to +1 for forward integration, and -1 for backward integration
%   stopAtRoot: stops integration when a root is encountered
%   stopAtR: stops integration when a particular value of r is encountered 
%   stopAfterR: stops integration after a particular value of r is encountered
%
%   Y is the last iteration of the log-derivative matrix.  rf is the final
%   separation at which integration stops.  nodes is the number of poles in
%   the log-derivative matrix that were encountered (only works for forward
%   integration).
%
%   If opt.output = true, then one can use [Y,rf,nodes,rout,Yout,Z,dbg] =
%   ... for additional output arguments.  rout is the vector of
%   separations, Yout is all of the log-derivative matrices, Z is the
%   transformation matrices used for constructing wavefunctions, and dbg is
%   a debugging structure.
%% Parse arguments
if nargin<4
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end

%% Prep
I = eye(ops.Nch);
if opt.direction>0
    ops = ops.int2spin;
    Ynew = real(sqrt(ops.potential(r(1),E))).*I;
elseif opt.direction<0
    r = flip(r);
    opt.blocks = rot90(numel(r)-opt.blocks+1,2);
    ops = ops.spin2int;
    Ynew = real(sqrt(ops.potential(r(1),E))).*I;
    Ynew = ops.U*Ynew*ops.U';
    ops = ops.rotate;
end
Ynew = sign(opt.direction)*diag(diag(Ynew));

[vecNew,eigNew] = eig(Ynew,'vector');
% eigNew = diag(eigNew);
eigNew = real(eigNew);
vecNew = gramschmidt(vecNew);
negEigNew = sum(eigNew<0);
eigAOld = eigNew;
eigANew = eigNew;
vecANew = vecNew;


if opt.output
    Y = zeros(ops.Nch,ops.Nch,numel(r));
    Y(:,:,1) = Ynew;
    Zout = zeros(ops.Nch,ops.Nch,numel(r));
    eigOut = zeros(ops.Nch,numel(r));
    eigOut(:,1) = eigNew;
end

breakFlag = false;
nodes = zeros(numel(r),1);
dbg.test = zeros(numel(r),1);
gcnt = 1;

for bb=1:size(opt.blocks,1)
    rb = r(opt.blocks(bb,1):opt.blocks(bb,2));
    h = (rb(2)-rb(1))/2;
%     Epot = repmat(E*I,1,1,numel(rb));
    M = ops.potential(rb,E);
    M2 = ops.potential(rb+h,E);
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
        [vecNew,eigNew] = eig(Ynew,'vector');
        eigNew = real(eigNew);
%         eigNew = diag(eigNew);
        vecNew = gramschmidt(vecNew);
        minEigNew = Inf;
        maxEigNew = 0;
        for kk=1:ops.Nch
            eigTest = abs(eigNew(kk));
            if eigTest<minEigNew
                minEigNew = eigTest;
            elseif eigTest>maxEigNew
                maxEigNew = eigTest;
            end
        end
        
        for knew=1:ops.Nch
            l2dist = 0;
            l2idx = 0;
            for kold=1:ops.Nch
                l2Test = abs(vecNew(:,knew)'*vecAOld(:,kold));
                if l2Test>=l2dist
                    l2dist = l2Test;
                    l2idx = kold;
                end
            end
            vecANew(:,l2idx) = vecNew(:,knew);
            eigANew(l2idx) = eigNew(knew);
        end
        
        if opt.output
            Y(:,:,gcnt+1) = Ynew;
            Zout(:,:,gcnt+1) = Zac*Zcb;
            eigOut(:,gcnt+1) = eigANew;
        end
        
        dnodes = 0;
        for kk=1:ops.Nch
%             if sign(eigAOld(kk)) ~= sign(eigANew(kk))
            if opt.direction>0 && eigAOld2(kk)<0 && eigAOld(kk)<0 && eigANew(kk)>0
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
    dbg.eigOut = eigOut(:,1:gcnt);
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
    
