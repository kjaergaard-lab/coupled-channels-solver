function Hdd=DipoleDipole(BV)
% DipoleDipole Generates dipole-dipole coupling matrix for coupled-channels calculations
%of akali-metal atom collisions from uncoupled basis vectors BV
Spin=0.5;
NumChannels=size(BV,1);
Hdd=zeros(NumChannels);
L=2;

for row=1:NumChannels,
    for col=1:NumChannels,        
        mI1f=BV(row,4);
        mI2f=BV(row,6);
        mI1i=BV(col,4);
        mI2i=BV(col,6);
        NSpinFactor=(mI1f==mI1i) && (mI2f==mI2i);
        if NSpinFactor==0,
            continue;
        end;
        
        Lf=BV(row,1);
        mLf=BV(row,2);
        Li=BV(col,1);
        mLi=BV(col,2);
        mS1f=BV(row,3);
        mS1i=BV(col,3);
        mS2f=BV(row,5);
        mS2i=BV(col,5);
             
        S1p=-1/sqrt(2)*sqrt(Spin*(Spin+1)-mS1i*(mS1i+1))*(mS1f==(mS1i+1));
        S1m=1/sqrt(2)*sqrt(Spin*(Spin+1)-mS1i*(mS1i-1))*(mS1f==(mS1i-1));       
        S1z=mS1i*(mS1i==mS1f);
           
        S2z=mS2i*(mS2i==mS2f);
        S2p=-1/sqrt(2)*sqrt(Spin*(Spin+1)-mS2i*(mS2i+1))*(mS2f==(mS2i+1));
        S2m=1/sqrt(2)*sqrt(Spin*(Spin+1)-mS2i*(mS2i-1))*(mS2f==(mS2i-1)); 
        
        Match=@(q) Lf>=abs(Li-L) && Lf<=(Li+L) && (mLi+q)==mLf;
        CG=@(q) sqrt((2*L+1)*(2*Li+1)./(4*pi*(2*Lf+1))).*ClebschGordan(Li,L,Lf,0,0,0)*ClebschGordan(Li,L,Lf,mLi,q,mLf);
%         CG=@(q) (-1).^(mLf+mLi).*sqrt((2*L+1)*(2*Lf+1)./(4*pi*(2*Li+1))).*ClebschGordan(Lf,L,Li,0,0,0)*ClebschGordan(Lf,L,Li,-mLf,q,-mLi);
        
        q=0;
        if Match(q),
            Ang0=CG(q);
        else
            Ang0=0;
        end;
        T0=Ang0.*(2*S1z*S2z+(S1p*S2m+S1m*S2p))/sqrt(6);
        
        q=1;
        if Match(q),
            Ang1p=CG(q);
        else
            Ang1p=0;
        end;
        T1p=-Ang1p.*(S1z*S2m+S1m*S2z)/sqrt(2);

        q=-1;
        if Match(q),
            Ang1m=CG(q);
        else
            Ang1m=0;
        end;
        T1m=-Ang1m.*(S1z*S2p+S1p*S2z)/sqrt(2);
        
        q=2;
        if Match(q),
            Ang2p=CG(q);
        else
            Ang2p=0;
        end;
        T2p=Ang2p.*(S1m*S2m);

        q=-2;
        if Match(q),
            Ang2m=CG(q);
        else
            Ang2m=0;
        end;
        T2m=Ang2m.*(S1p*S2p);
        
        
        Hdd(row,col)=sqrt(4*pi/5)*(T0+T1p+T1m+T2p+T2m).*NSpinFactor;
    end;
end;
        
        