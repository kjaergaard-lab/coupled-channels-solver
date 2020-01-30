function U=RbBOPotentialsv2(r,sel)
% v2 update uses retardation values from PRA 50(4) 3096
% Energies in cm^-1, lengths in angstroms
if sel==0,
    Rinner=3.126;
    Router=11;

    %% Inner region
    A=-0.638904880e4;
    B=0.112005361e7;
    Ns=4.53389;

    Uinner=A+B./r.^Ns;


    %% Intermediate region
    b=-0.13;
    Rm=4.209912760;
    a(1)=-3993.592873;
    a(2)=0;
    a(3)=0.282069372972346137e5;
    a(4)=0.560425000209256905e4;
    a(5)=-0.423962138510562945e5;
    a(6)=-0.598558066508841584e5;
    a(7)=-0.162613532034769596e5;
    a(8)=-0.405142102246254944e5;
    a(9)=0.195237415352729586e6;
    a(10)=0.413823663033582852e6;
    a(11)=-0.425543284828921501e7;
    a(12)=0.546674790157210198e6;
    a(13)=0.663194778861331940e8;
    a(14)=-0.558341849704095051e8;
    a(15)=-0.573987344918535471e9;
    a(16)=0.102010964189156187e10;
    a(17)=0.300040150506311035e10;
    a(18)=-0.893187252759830856e10;
    a(19)=-0.736002541483347511e10;
    a(20)=0.423130460980355225e11;
    a(21)=-0.786351477693491840e10;
    a(22)=-0.102470557344862152e12;
    a(23)=0.895155811349267578e11;
    a(24)=0.830355322355692902e11;
    a(25)=-0.150102297761234375e12;
    a(26)=0.586778574293387070e11;

    x=(r-Rm)./(r+b.*Rm);
    Nterms=numel(a);
    Uir=zeros(size(r));
    for nn=1:Nterms,
        Uir=Uir+a(nn).*x.^(nn-1);
    end;

    %% Outer Region
    Uinf=0;
    C6=0.2270032e8;
    C8=0.7782886e9;
    C10=0.2868869e11;
    C26=0.2819810e26;
    Aex=0.1317786e5;
    Gamma=5.317689;
    Beta=2.093816;
    
    rtmp=r(:);
    rtmp=rtmp(rtmp>Router);
    [f6,f8,f10]=RbRetardation(rtmp);
    f6=reshape([zeros(numel(r)-numel(rtmp),1);f6],size(r));
    f8=reshape([zeros(numel(r)-numel(rtmp),1);f8],size(r));
    f10=reshape([zeros(numel(r)-numel(rtmp),1);f10],size(r));

%     Uouter=Uinf-C6./r.^6-C8./r.^8-C10./r.^10-C26./r.^26-Aex.*r.^Gamma.*exp(-Beta.*r);
    Uouter=Uinf-C6.*f6./r.^6-C8.*f8./r.^8-C10.*f10./r.^10-C26./r.^26-Aex.*r.^Gamma.*exp(-Beta.*r);
else
    Rinner=5.07;
    Router=11;

    %% Inner region
    Ns=4.5338950;
    A=-0.619088543e3;
    B=0.956231677e6;

    Uinner=A+B./r.^Ns;


    %% Intermediate region
    b=-0.33;
    Rm=6.0933451;
    a(1)=-241.503352;
    a(2)=-0.672503402304666542;
    a(3)=0.195494577140503543e4;
    a(4)=-0.141544168453406223e4;
    a(5)=-0.221166468149940465e4;
    a(6)=0.165443726445793004e4;
    a(7)=-0.596412188910614259e4;
    a(8)=0.654481694231538040e4;
    a(9)=0.261413416681972012e5;
    a(10)=-0.349701859112702878e5;
    a(11)=-0.328185277155018630e5;
    a(12)=0.790208849885562522e5;
    a(13)=-0.398783520249289213e5;
    
    x=(r-Rm)./(r+b.*Rm);
    Nterms=numel(a);
    Uir=zeros(size(r));
    for nn=1:Nterms,
        Uir=Uir+a(nn).*x.^(nn-1);
    end;

    %% Outer Region
    Uinf=0;
    C6=0.2270032e8;
    C8=0.7782886e9;
    C10=0.2868869e11;
    C26=0.2819810e26;
    Aex=0.1317786e5;
    Gamma=5.317689;
    Beta=2.093816;
    
    rtmp=r(:);
    rtmp=rtmp(rtmp>Router);
    [f6,f8,f10]=RbRetardation(rtmp);
    f6=reshape([ones(numel(r)-numel(rtmp),1);f6],size(r));
    f8=reshape([ones(numel(r)-numel(rtmp),1);f8],size(r));
    f10=reshape([ones(numel(r)-numel(rtmp),1);f10],size(r));

%     Uouter=Uinf-C6./r.^6-C8./r.^8-C10./r.^10-C26./r.^26+Aex.*r.^Gamma.*exp(-Beta.*r);
    Uouter=Uinf-C6.*f6./r.^6-C8.*f8./r.^8-C10.*f10./r.^10-C26./r.^26+Aex.*r.^Gamma.*exp(-Beta.*r);
end;

%% All together
U=Uinner.*(r<Rinner) + Uir.*((r>=Rinner) & (r<=Router)) + Uouter.*(r>Router);


end

