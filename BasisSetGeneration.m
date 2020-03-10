function BasisSetGeneration(AtomPair,LVec,filename)
% BasisSetGeneration Creates BasisSetFile for use in MultiChannel
%   BasisSetGeneration(AtomPair,LVec,filename)
%   AtomPair is a string and one of 'RbRb', 'KK', or 'KRb'.
%   LVec defines the range of angular momentum quantum numbers to include,
%   as for instance [0:2], or [1,3,5]
%   filename is the name of the file that is saved.  Default is [AtomPair,
%   ' Basis Set L=LVec']
%   See also MultiChannel, ClebschGordan

%% Preamble
constants;
if nargin==2
    Lstr=ProcessLstr(LVec);
    filename=[AtomPair,' Basis Set ',Lstr];
end

%% Non-user settable parameters
switch AtomPair
    case 'KRb'
        PotentialFunc=@KRbPotentialDD;  %Function describing the interaction potential

        %Angular momentum settings
        Spin=0.5;   %Electron spin
        NSpin1=4;   %Spin of first atom
        NSpin2=1.5; %Spin of second atom

        %Set constant values
        IdenticalParticles=0;   %Important parameter for identical particles.  +1 for bosons, -1 for fermions, 0 for distinct particles
        mass=mK*mRb./(mK+mRb);  %Reduced mass
        Ahfs1=Ahfs_K;   %Hyperfine constant for atom 1
        Ahfs2=Ahfs_Rb;  %Hyperfine constant for atom 2
        gS_1=gS_K;
        gS_2=gS_Rb;
        gI_1=gI_K;
        gI_2=gI_Rb;

    case 'RbRb'
        PotentialFunc=@RbRbPotentialDD;

        %Angular momentum settings
        Spin=0.5;   %Electron spin
        NSpin1=1.5;   %Spin of first atom
        NSpin2=1.5; %Spin of second atom

        %Set constant values
        IdenticalParticles=1;   %Important parameter for identical particles.  +1 for bosons, -1 for fermions, 0 for distinct particles
        mass=mRb/2;
        Ahfs1=Ahfs_Rb;
        Ahfs2=Ahfs_Rb;
        gS_1=gS_Rb;
        gS_2=gS_Rb;
        gI_1=gI_Rb;
        gI_2=gI_Rb;
        
    case 'KK'
        PotentialFunc=@KKPotentialDD;

        %Angular momentum settings
        Spin=0.5;   %Electron spin
        NSpin1=4;   %Spin of first atom
        NSpin2=4; %Spin of second atom

        %Set constant values
        IdenticalParticles=-1;   %Important parameter for identical particles.  +1 for bosons, -1 for fermions, 0 for distinct particles
        mass=mK/2;
        Ahfs1=Ahfs_K;
        Ahfs2=Ahfs_K;
        gS_1=gS_K;
        gS_2=gS_K;
        gI_1=gI_K;
        gI_2=gI_K;  
      
    case '41K41K'
        PotentialFunc=@KKPotentialDD;

        %Angular momentum settings
        Spin=0.5;   %Electron spin
        NSpin1=1.5;   %Spin of first atom
        NSpin2=1.5; %Spin of second atom

        %Set constant values
        IdenticalParticles=1;   %Important parameter for identical particles.  +1 for bosons, -1 for fermions, 0 for distinct particles
        mass=mK41/2;
        Ahfs1=Ahfs_K41;
        Ahfs2=Ahfs_K41;
        gS_1=gS_K41;
        gS_2=gS_K41;
        gI_1=gI_K41;
        gI_2=gI_K41;

    case '41K87Rb'
        PotentialFunc=@KRbPotentialDD;  %Function describing the interaction potential

        %Angular momentum settings
        Spin=0.5;   %Electron spin
        NSpin1=1.5;   %Spin of first atom
        NSpin2=1.5; %Spin of second atom

        %Set constant values
        IdenticalParticles=0;   %Important parameter for identical particles.  +1 for bosons, -1 for fermions, 0 for distinct particles
        mass=mK41*mRb./(mK41+mRb);  %Reduced mass
        Ahfs1=Ahfs_K41;   %Hyperfine constant for atom 1
        Ahfs2=Ahfs_Rb;  %Hyperfine constant for atom 2
        gS_1=gS_K41;
        gS_2=gS_Rb;
        gI_1=gI_K41;
        gI_2=gI_Rb;
        
    case '85Rb85Rb'
        PotentialFunc=@RbRbPotentialDD;

        %Angular momentum settings
        Spin=0.5;   %Electron spin
        NSpin1=2.5;   %Spin of first atom
        NSpin2=2.5; %Spin of second atom

        %Set constant values
        IdenticalParticles=1;   %Important parameter for identical particles.  +1 for bosons, -1 for fermions, 0 for distinct particles
        mass=mRb85/2;
        Ahfs1=Ahfs_Rb85;
        Ahfs2=Ahfs_Rb85;
        gS_1=gS_Rb85;
        gS_2=gS_Rb85;
        gI_1=gI_Rb85;
        gI_2=gI_Rb85;
        
    case '39K39K'
        PotentialFunc=@KKPotentialDD;

        %Angular momentum settings
        Spin=0.5;   %Electron spin
        NSpin1=1.5;   %Spin of first atom
        NSpin2=1.5; %Spin of second atom

        %Set constant values
        IdenticalParticles=1;   %Important parameter for identical particles.  +1 for bosons, -1 for fermions, 0 for distinct particles
        mass=mK39/2;
        Ahfs1=Ahfs_K39;
        Ahfs2=Ahfs_K39;
        gS_1=gS_K39;
        gS_2=gS_K39;
        gI_1=gI_K39;
        gI_2=gI_K39;
        
    otherwise
        error('Unsupported atom pair!  No file generated.');
end

NumSpinStates=(2*Spin+1);
NumNSpinStates2=(2*NSpin2+1);
NumNSpinStates1=(2*NSpin1+1);
NumTSpinStates2=(NumSpinStates)*(NumNSpinStates2);
NumTSpinStates1=(NumSpinStates)*(NumNSpinStates1);
NumTSpinStates=NumTSpinStates2*NumTSpinStates1;
NumLStates=sum(2*LVec+1);

NumChannels=NumLStates*NumTSpinStates;
BV1=zeros(NumChannels,6);
BV2=BV1;
BV3=BV1;
BV4=BV1;
BV5=BV1;
BVint=zeros(NumChannels,4);
I=eye(NumChannels);


%% Basis state definitions
% Fully uncoupled basis vectors BV1
% |L mL mSK mIK mSRb mIRb> = |L mL mS1 mI1 mS2 mI2>
count=1;
for n1=1:numel(LVec)
    for n2=1:(2*LVec(n1)+1)
        for n3=1:NumSpinStates
            for n4=1:NumNSpinStates1
                for n5=1:NumSpinStates
                    for n6=1:NumNSpinStates2
                        mL=-LVec(n1)+(n2-1);
                        mS1=-Spin+(n3-1);
                        mI1=-NSpin1+(n4-1);
                        mS2=-Spin+(n5-1);
                        mI2=-NSpin2+(n6-1);
                        BV1(count,:)=[LVec(n1),mL,mS1,mI1,mS2,mI2];
                        count=count+1;
                    end
                end
            end
        end
    end
end
LMat=diag(BV1(:,1));

% Spin coupled basis vectors BV2 S=S1+S2
% |L mL S mS mIK mIRb> = |L mL S mS mI1 mI2>
count=1;
for n1=1:numel(LVec)
    for n2=1:(2*LVec(n1)+1)
        SVec=0:1;
        for n3=1:numel(SVec)
            for n4=1:(2*SVec(n3)+1)
                for n5=1:NumNSpinStates1
                    for n6=1:NumNSpinStates2
                        mL=-LVec(n1)+(n2-1);
                        S=SVec(n3);
                        mS=-SVec(n3)+(n4-1);
                        mI1=-NSpin1+(n5-1);
                        mI2=-NSpin2+(n6-1);
                        BV2(count,:)=[LVec(n1),mL,S,mS,mI1,mI2];
                        count=count+1;
                    end
                end
            end
        end
    end
end

% Spin-nucleus coupled basis vectors BV3 F1=S1+I1 and F2=S2+I2
% |L mL FK mFK FRb mFRb> = |L mL F1 mF1 F2 mF2>
count=1;
for n1=1:numel(LVec)
    for n2=1:(2*LVec(n1)+1)
        FVec1=abs(Spin-NSpin1):(Spin+NSpin1);
        for n3=1:numel(FVec1)
            for n4=1:(2*FVec1(n3)+1)
                FVec2=abs(Spin-NSpin2):(Spin+NSpin2);
                for n5=1:numel(FVec2)
                    for n6=1:(2*FVec2(n5)+1)
                        mL=-LVec(n1)+(n2-1);
                        F1=FVec1(n3);
                        mF1=-FVec1(n3)+(n4-1);
                        F2=FVec2(n5);
                        mF2=-FVec2(n5)+(n6-1);
                        BV3(count,:)=[LVec(n1),mL,F1,mF1,F2,mF2];
                        count=count+1;
                    end
                end
            end
        end
    end
end

% Coupled total spin vectors BV4 F1+F2=F
% |L mL F1 F2 F mF>
count=1;
for n1=1:numel(LVec)
    for n2=1:(2*LVec(n1)+1)
         FVec1=abs(Spin-NSpin1):(Spin+NSpin1);
         for n3=1:numel(FVec1)
            FVec2=abs(Spin-NSpin2):(Spin+NSpin2);
            for n4=1:numel(FVec2)
                FVec=abs(FVec1(n3)-FVec2(n4)):(FVec1(n3)+FVec2(n4));
                for n5=1:numel(FVec)
                    for n6=1:(2*FVec(n5)+1)
                        L=LVec(n1);
                        mL=-L+(n2-1);
                        F1=FVec1(n3);
                        F2=FVec2(n4);
                        F=FVec(n5);
                        mF=-F+(n6-1);
                        BV4(count,:)=[L,mL,F1,F2,F,mF];
                        count=count+1;
                    end
                end
            end
         end
    end
end

% Fully coupled basis vectors BV5 J=L+F=L+F1+F2=L+S1+S2+I1+I2
% |L F1 F2 F J mJ>
count=1;
for n1=1:numel(LVec)
    FVec1=abs(Spin-NSpin1):(Spin+NSpin1);
    for n2=1:numel(FVec1)
        FVec2=abs(Spin-NSpin2):(Spin+NSpin2);
        for n3=1:numel(FVec2)
            FVec=abs(FVec1(n2)-FVec2(n3)):(FVec1(n2)+FVec2(n3));
            for n4=1:numel(FVec)
                JVec=abs(FVec(n4)-LVec(n1)):(FVec(n4)+LVec(n1));
                for n5=1:numel(JVec)
                    for n6=1:(2*JVec(n5)+1)
                        L=LVec(n1);
                        F1=FVec1(n2);
                        F2=FVec2(n3);
                        F=FVec(n4);
                        J=JVec(n5);
                        mJ=-JVec(n5)+(n6-1);
                        BV5(count,:)=[L,F1,F2,F,J,mJ];
                        count=count+1;
                    end
                end
            end
        end
    end
end

% Basis of internal states, labelled from lowest energy BVint
% |L mL Int1 Int2>
count=1;
for n1=1:numel(LVec)
    for n2=1:(2*LVec(n1)+1)
        for n3=1:(NumTSpinStates1)
            for n4=1:(NumTSpinStates2)
                    mL=-LVec(n1)+(n2-1);
                    BVint(count,:)=[LVec(n1),mL,n3,n4];
                    count=count+1;
            end
        end
    end
end

%% Assign mF values to internal states such that internal states are labelled in increasing energy
FSpin1=NSpin1+[Spin,-Spin];
FSpin2=NSpin2+[Spin,-Spin];
if Ahfs1<0
%     IntMValues(:,1)=(BVint(:,3)<=(2*FSpin1(1)+1)).*(BVint(:,3)-(FSpin1(1)+1))+(BVint(:,3)>(2*FSpin1(1)+1)).*(2*FSpin1(1)+2+FSpin1(2)-BVint(:,3));
    IntMValues(:,1)=(BVint(:,3)<=(2*FSpin1(1)+1)).*(BVint(:,3)-(FSpin1(1)+1))+(BVint(:,3)>(2*FSpin1(1)+1)).*((2*FSpin1(1)+1)+(FSpin1(2)+1)-BVint(:,3));
else
    IntMValues(:,1)=(BVint(:,3)<=(2*FSpin1(2)+1)).*((FSpin1(2)+1)-BVint(:,3))+(BVint(:,3)>(2*FSpin1(2)+1)).*(BVint(:,3)-(2*FSpin1(2)+1)-(FSpin1(1)+1));
end

if Ahfs2<0
    IntMValues(:,2)=(BVint(:,4)<=(2*FSpin2(1)+1)).*(BVint(:,4)-(FSpin2(1)+1))+(BVint(:,4)>(2*FSpin2(1)+1)).*((2*FSpin2(1)+1)+(FSpin2(2)+1)-BVint(:,4));
else
    IntMValues(:,2)=(BVint(:,4)<=(2*FSpin2(2)+1)).*((FSpin2(2)+1)-BVint(:,4))+(BVint(:,4)>(2*FSpin2(2)+1)).*(BVint(:,4)-(2*FSpin2(2)+1)-(FSpin2(1)+1));
end


%% Basis transformations

% Transformation matrix from fully uncoupled to spin-coupled basis
BT21=zeros(NumChannels);
for row=1:NumChannels
    for col=1:NumChannels
        NSpin1Match=BV1(col,4)==BV2(row,5);
        NSpin2Match=BV1(col,6)==BV2(row,6);
        SpinSumMatch=(BV1(col,3)+BV1(col,5))==BV2(row,4);
        mLMatch=BV1(col,2)==BV2(row,2);
        LMatch=BV1(col,1)==BV2(row,1);
        if NSpin1Match && NSpin2Match && SpinSumMatch && mLMatch && LMatch
            BT21(row,col)=ClebschGordan(Spin,Spin,BV2(row,3),BV1(col,3),BV1(col,5),BV2(row,4));
        end
    end
end

% Transformation matrix from fully uncoupled to spin-nucleus coupled basis
BT31=zeros(NumChannels);
for row=1:NumChannels
    for col=1:NumChannels
        SpinSumMatch1=(BV1(col,3)+BV1(col,4))==BV3(row,4);
        SpinSumMatch2=(BV1(col,5)+BV1(col,6))==BV3(row,6);
        mLMatch=BV1(col,2)==BV3(row,2);
        LMatch=BV1(col,1)==BV3(row,1);
        if SpinSumMatch1 && SpinSumMatch2 && mLMatch && LMatch,
            tmp1=ClebschGordan(Spin,NSpin1,BV3(row,3),BV1(col,3),BV1(col,4),BV3(row,4));
            tmp2=ClebschGordan(Spin,NSpin2,BV3(row,5),BV1(col,5),BV1(col,6),BV3(row,6));
            BT31(row,col)=tmp1.*tmp2;
        end
    end
end

% Transformation matrix from spin-coupled basis to spin-nucleus coupled
% basis
BT32=BT31*(BT21');

% Transformation matrix from spin-nucleus basis to total spin coupled
BT43=zeros(NumChannels);
for row=1:NumChannels
    for col=1:NumChannels
        mLMatch=BV3(col,2)==BV4(row,2);
        LMatch=BV3(col,1)==BV4(row,1);
        FMatch1=BV3(col,3)==BV4(row,3);
        FMatch2=BV3(col,5)==BV4(row,4);
        FMatch=FMatch1 && FMatch2;
        mFSumMatch=(BV3(col,4)+BV3(col,6))==BV4(row,6);
        if mLMatch && LMatch && mFSumMatch && FMatch
            BT43(row,col)=ClebschGordan(BV3(col,3),BV3(col,5),BV4(row,5),BV3(col,4),BV3(col,6),BV4(row,6));
        end
    end
end

% Transformation matrix from total-spin coupled basis to fully coupled
BT54=zeros(NumChannels);
for row=1:NumChannels
    for col=1:NumChannels
        LMatch=BV4(col,1)==BV5(row,1);
        FMatch1=BV4(col,3)==BV5(row,2);
        FMatch2=BV4(col,4)==BV5(row,3);        
        FTotalMatch=BV4(col,5)==BV5(row,4);
        FMatch=FMatch1 && FMatch2 && FTotalMatch;
        mJSumMatch=(BV4(col,2)+BV4(col,6))==BV5(row,6);
        if LMatch && FMatch && mJSumMatch
            BT54(row,col)=ClebschGordan(BV4(col,1),BV4(col,5),BV5(row,5),BV4(col,2),BV4(col,6),BV5(row,6));
        end
    end
end

% Transformation matrix from electron spin-coupled basis to total spin coupled
BT42=BT43*BT32;

% save(filename,'LMat','BV1','BV2','BV3','BV4','BV5','BT21','BT31','BT32','BT43','BT54','BT42');
save(filename);

end

function Lstr=ProcessLstr(LVec)
% Produces a string from LVec that compresses concurrent numbers
% So LVec=[0:2,5,8:10] becomes 'L=0-2,5,8-10'

LVec=sort(unique(LVec));
if numel(LVec)==2
    if diff(LVec)==1
        Lstr=['L=',num2str(LVec(1)),'-',num2str(LVec(2))];
    else
        Lstr=['L=',num2str(LVec(1)),',',num2str(LVec(2))];
    end
else
    Lstr=['L=',num2str(LVec(1))];
end

for nn=2:numel(LVec)-1
    dL=diff(LVec(nn-1:nn+1));
    if dL(1)==1 && dL(2)~=1
        Lstr=[Lstr,'-',num2str(LVec(nn))];
    elseif dL(1)~=1 && dL(2)==1
        Lstr=[Lstr,',',num2str(LVec(nn))];
    elseif dL(1)~=1 && dL(2)~=1
        Lstr=[Lstr,',',num2str(LVec(nn))];
    end
    
    if nn==numel(LVec)-1
        if dL(2)==1
           Lstr=[Lstr,'-',num2str(LVec(nn+1))];
        else
            Lstr=[Lstr,',',num2str(LVec(nn+1))];
        end;           
    end
end
end