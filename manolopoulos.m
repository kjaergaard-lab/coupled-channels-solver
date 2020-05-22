function mout = manolopoulos(E,ops,basis,opt)


%% Prep
I = eye(ops.Nch);
Ynew = sqrt(ops.potential(opt.rmin,E)).*I;
if opt.getwf
    Zout = zeros(ops.Nch);
    rout = opt.rmin;
end

%% Calculate adiabatic potentials
ra = linspace(opt.rmin,opt.rmax,5e3);
f = ops.potential(ra,E);
adiabat = zeros(ops.Nch,numel(ra));
for kk=1:numel(ra)
    adiabat(:,kk) = sort(eig(f(:,:,kk)));
end

changeFlag = false;
numBlocks = floor(opt.rmax/opt.blocksize);

for jj=1:numBlocks
    rstart = opt.rmin+(jj-1)*opt.blocksize;
    rend = rstart+opt.blocksize;
    
    kmax = adiabat(1,ra>=rstart & ra<=rend);
    kmax = kmax(kmax<0);
    kmax = max(sqrt(abs(kmax)));
    if ~isempty(kmax)
        dr = opt.drscale*2*pi./kmax;
        dr = max(min(dr,opt.drmax),opt.drmin);
        dr = min(dr,opt.blocksize/2);
    end
    
    N = max(ceil((rend-rstart)/dr),2);
    rb = linspace(rstart,rend,N);
    h = (rb(2)-rb(1))/2;
    
    if opt.getwf
        Z = zeros(ops.Nch,ops.Nch,N);
    end

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
        
        if opt.getwf
            Z(:,:,nn+1) = Zac*Zcb;
            rout = [rout rb(2:end)]; %#ok<AGROW>
        end
    end
    
    if opt.getwf
        Zout(:,:,end+1:end+size(Z,3)-1) = Z(:,:,2:end);
    end
    
    if ~changeFlag && rend>=10
        changeFlag = true;
        ops = ops.rotate;
        Ynew = ops.U*Ynew*ops.U';
        if nargout>1
            for kk=1:size(Zout,3)
                Zout(:,:,kk) = ops.U*Zout(:,:,kk)*ops.U';
            end
        end
    end
end

%% Calculate S-matrix
b = rb(end);
Eabs = diag((E+ops.refE)*I-ops.Hint);
idx = Eabs>0;             %Determine which channels are open
Eabs = diag(Eabs);
k = sqrt(Eabs(idx,idx));
L = ops.L(idx,idx);
Ynew = Ynew(idx,idx);
Stmp = (k*RiccatiBessel(b,k,L,2,1)-RiccatiBessel(b,k,L,2,0)*Ynew)/(k*RiccatiBessel(b,k,L,1,1)-RiccatiBessel(b,k,L,1,0)*Ynew);
S = I;
S(idx,idx) = Stmp;

%% Calculate wavefunctions
if opt.getwf
    u = zeros(ops.Nch,size(Zout,3));
    if basis.symmetry == 0
        U = I;
    else
        U = 0*I;
        for row=1:size(ops.Nch,1)
            for col=1:size(ops.Nch,2)
                if all(basis.bvint(row,1:2) == basis.bvint(col,1:2),2)
                    if all(basis.bvint(row,5:6) == basis.bvint(col,5:6),2)
                        U(row,col) = U(row,col) + (2*(1+(basis.bvint(row,5)==basis.bvint(row,6))))^-0.5;
                    end

                    if all(basis.bvint(row,5:6) == basis.bvint(col,[6,5]),2)
                        U(row,col) = U(row,col) + (-1)^basis.bvint(row,1)*(2*(1+(basis.bvint(row,5)==basis.bvint(row,6))))^-0.5;
                    end
                end
            end
        end
    end
    
    phi = I;
    phi(idx,idx) = RiccatiBessel(b,k,L,2,0)-RiccatiBessel(b,k,L,1,0)*Stmp;
    C = U*diag(-sqrt(pi*(2*L+I)/k)*(1i)^(L-I));
    u(:,end) = phi*C;
    for nn=size(Zout,3):-1:2
        phi = Zout(:,:,nn)*phi;
        u(:,nn-1) = phi*C;
    end
    wf.r = rout;
    wf.u = u;
end

if opt.getwf
    mout.S = S;
    mout.wf = wf;
else
    mout = S;
end

end
    
