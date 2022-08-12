classdef atompairbasis < matlab.mixin.Copyable
    %ATOMPAIRBASIS Provides definitions of various basis vectors useful for
    %calculating scattering and bound state properties
    properties
        atoms       %A 1x2 cell array of ATOMPROPERTIES objects
        pair        %A character string indicating what BO potential to use
        
        espin       %The electron spin -- always 0.5
        nspin       %A 1x2 array of nuclear spins
        
        symmetry    %The symmetry of the atom pair: 0 for distinguishable particles, and +/-1 for bosons/fermions
        mass        %The reduced mass of the atom pair
        Ahfs        %A 1x2 array of hyperfine constants in MHz
        gS          %A 1x2 array of electron g-factors
        gI          %A 1x2 array of nuclear g-factors
        
        scale       %A scale factor converting cm^{-1} to inverse square wavenumbers in Angstrom^{-2}
        
        Nchannels   %The number of channels
        L           %A diagonal matrix of L values
        bv1         %Basis vectors in uncoupled basis: [L mL mS1 mI1 mS2 mI2]
        bv2         %Basis vectors in electron spin coupled basis S = S1+S2: [L mL S mS mI1 mI2]
        bv3         %Basis vectors in electron+nuclear spin coupled basis F = S+I: [L mL F1 mF1 F2 mF2]
        bvint       %Basis vectors in the internal basis: [L mL mF1 mF2 int1 int2]
        
        bt21        %Transformation matrix from basis 1 to basis 2
        bt31        %Transformation matrix from basis 1 to basis 3
        bt32        %Transformation matrix from basis 2 to basis 3
    end
    
    methods
        function self = atompairbasis(pair,Lvec)
            %ATOMPAIRBASIS Constructs an instance of the ATOMPAIRBASIS
            %class
            %
            %   basis = ATOMPAIRBASIS(pair,Lvec) constructs an object for
            %   the given pair (a character string) and vector of L values
            self.espin = 0.5;
            switch pair
                case {'KRb','K40Rb87'}
                    self.pair = 'krb';
                    atoms = {scattconst.K40,scattconst.Rb87};
                    self.symmetry = 0;
                case {'RbRb','Rb87Rb87'}
                    self.pair = 'rbrb';
                    atoms = {scattconst.Rb87,scattconst.Rb87};
                    self.symmetry = 1;
                case {'KK','K40K40'}
                    self.pair = 'kk';
                    atoms = {scattconst.K40,scattconst.K40};
                    self.symmetry = -1;
                case {'K41K41'}
                    self.pair = 'kk';
                    atoms = {scattconst.K41,scattconst.K41};
                    self.symmetry = 1;
                case {'K40K41'}
                    self.pair = 'kk';
                    atoms = {scattconst.K40,scattconst.K41};
                    self.symmetry = 0;
                case {'K39K39'}
                    self.pair = 'kk';
                    atoms = {scattconst.K39,scattconst.K39};
                    self.symmetry = 1;
                case {'K39K41'}
                    self.pair = 'kk';
                    atoms = {scattconst.K39,scattconst.K41};
                    self.symmetry = 0;
                case {'K39K40'}
                    self.pair = 'kk';
                    atoms = {scattconst.K39,scattconst.K40};
                    self.symmetry = 0;
                case {'K41Rb87'}
                    self.pair = 'krb';
                    atoms = {scattconst.K41,scattconst.Rb87};
                    self.symmetry = 0;
                case {'K39Rb87'}
                    self.pair = 'krb';
                    atoms = {scattconst.K39,scattconst.Rb87};
                    self.symmetry = 0;
                case {'K39Rb85'}
                    self.pair = 'krb';
                    atoms = {scattconst.K39,scattconst.Rb85};
                    self.symmetry = 0;
                case {'K40Rb85'}
                    self.pair = 'krb';
                    atoms = {scattconst.K40,scattconst.Rb85};
                    self.symmetry = 0;
                case {'K41Rb85'}
                    self.pair = 'krb';
                    atoms = {scattconst.K41,scattconst.Rb85};
                    self.symmetry = 0;
                case {'Rb85Rb85'}
                    self.pair = 'rbrb';
                    atoms = {scattconst.Rb85,scattconst.Rb85};
                    self.symmetry = 1;
                case {'Rb85Rb87'}
                    self.pair = 'rbrb';
                    atoms = {scattconst.Rb85,scattconst.Rb87};
                    self.symmetry = 0;
                otherwise
                    error('Atom pair not supported');
            end
            self.nspin = [atoms{1}.nspin,atoms{2}.nspin];
            self.mass = atoms{1}.mass*atoms{2}.mass/(atoms{1}.mass+atoms{2}.mass);
            self.scale = scattconst.cm2K*scattconst.K2A(self.mass);
            self.Ahfs = [atoms{1}.Ahfs,atoms{2}.Ahfs];
            self.gS = [atoms{1}.gS,atoms{2}.gS];
            self.gI = [atoms{1}.gI,atoms{2}.gI];
            self.atoms = atoms;
            self.Nchannels = sum(2*Lvec+1)*(2*self.espin+1)^2*prod(2*self.nspin+1);
            
            self.makeBasisVectors(Lvec);
            self.makeBasisTransformations;
            
        end
        
        function self = makeBasisVectors(self,LVec)
            %MAKEBASISVECTORS Creates basis vectors given a vector of L
            %values
            %
            %   self = MAKEBASISVECTORS(LVec) creates basis vectors for a
            %   given vector of L values
            self.bv1 = zeros(self.Nchannels,6);
            self.bv2 = zeros(self.Nchannels,6);
            self.bv3 = zeros(self.Nchannels,6);
            self.bvint = zeros(self.Nchannels,6);
            %% Basis state definitions
            % Fully uncoupled basis vectors BV1
            % |L mL mSK mIK mSRb mIRb> = |L mL mS1 mI1 mS2 mI2>
            count = 1;
            for n1=1:numel(LVec)
                for n2=1:(2*LVec(n1)+1)
                    for n3=1:(2*self.espin+1)
                        for n4=1:(2*self.nspin(1)+1)
                            for n5=1:(2*self.espin+1)
                                for n6=1:(2*self.nspin(2)+1)
                                    mL = -LVec(n1)+(n2-1);
                                    mS1 = -self.espin+(n3-1);
                                    mI1 = -self.nspin(1)+(n4-1);
                                    mS2 = -self.espin+(n5-1);
                                    mI2 = -self.nspin(2)+(n6-1);
                                    self.bv1(count,:) = [LVec(n1),mL,mS1,mI1,mS2,mI2];
                                    count = count+1;
                                end
                            end
                        end
                    end
                end
            end
            self.L = diag(self.bv1(:,1));
            
            % Spin coupled basis vectors BV2 S=S1+S2
            % |L mL S mS mIK mIRb> = |L mL S mS mI1 mI2>
            count = 1;
            for n1=1:numel(LVec)
                for n2=1:(2*LVec(n1)+1)
                    SVec = 0:1;
                    for n3=1:numel(SVec)
                        for n4=1:(2*SVec(n3)+1)
                            for n5=1:(2*self.nspin(1)+1)
                                for n6=1:(2*self.nspin(2)+1)
                                    mL = -LVec(n1)+(n2-1);
                                    S = SVec(n3);
                                    mS = -SVec(n3)+(n4-1);
                                    mI1 = -self.nspin(1)+(n5-1);
                                    mI2 = -self.nspin(2)+(n6-1);
                                    self.bv2(count,:) = [LVec(n1),mL,S,mS,mI1,mI2];
                                    count = count+1;
                                end
                            end
                        end
                    end
                end
            end
            
            % Spin-nucleus coupled basis vectors BV3 F1=S1+I1 and F2=S2+I2
            % |L mL FK mFK FRb mFRb> = |L mL F1 mF1 F2 mF2>
            count = 1;
            for n1=1:numel(LVec)
                for n2=1:(2*LVec(n1)+1)
                    FVec1 = abs(self.espin-self.nspin(1)):(self.espin+self.nspin(1));
                    for n3=1:numel(FVec1)
                        for n4=1:(2*FVec1(n3)+1)
                            FVec2 = abs(self.espin-self.nspin(2)):(self.espin+self.nspin(2));
                            for n5=1:numel(FVec2)
                                for n6=1:(2*FVec2(n5)+1)
                                    mL = -LVec(n1)+(n2-1);
                                    F1 = FVec1(n3);
                                    mF1 = -FVec1(n3)+(n4-1);
                                    F2 = FVec2(n5);
                                    mF2 = -FVec2(n5)+(n6-1);
                                    self.bv3(count,:) = [LVec(n1),mL,F1,mF1,F2,mF2];
                                    count = count+1;
                                end
                            end
                        end
                    end
                end
            end
            
            % Basis of internal states, labelled from lowest energy BVint
            % |L mL Int1 Int2>
            count = 1;
            FSpin1 = self.nspin(1)+self.espin*[1,-1];
            FSpin2 = self.nspin(2)+self.espin*[1,-1];
            for n1=1:numel(LVec)
                for n2=1:(2*LVec(n1)+1)
                    for n3=1:((2*self.espin+1)*(2*self.nspin(1)+1))
                        for n4=1:((2*self.espin+1)*(2*self.nspin(2)+1))
                            mL = -LVec(n1)+(n2-1);
                            % Assign mF values to internal states such that internal states are labelled in increasing energy
                            if self.Ahfs(1)<0
                                mF3 = (n3<=(2*FSpin1(1)+1)).*(n3-(FSpin1(1)+1))+(n3>(2*FSpin1(1)+1)).*((2*FSpin1(1)+1)+(FSpin1(2)+1)-n3);
                            else
                                mF3 = (n3<=(2*FSpin1(2)+1)).*((FSpin1(2)+1)-n3)+(n3>(2*FSpin1(2)+1)).*(n3-(2*FSpin1(2)+1)-(FSpin1(1)+1));
                            end
                            
                            if self.Ahfs(2)<0
                                mF4 = (n4<=(2*FSpin2(1)+1)).*(n4-(FSpin2(1)+1))+(n4>(2*FSpin2(1)+1)).*((2*FSpin2(1)+1)+(FSpin2(2)+1)-n4);
                            else
                                mF4 = (n4<=(2*FSpin2(2)+1)).*((FSpin2(2)+1)-n4)+(n4>(2*FSpin2(2)+1)).*(n4-(2*FSpin2(2)+1)-(FSpin2(1)+1));
                            end
                            
                            self.bvint(count,:) = [LVec(n1),mL,mF3,mF4,n3,n4];
                            count = count+1;
                        end
                    end
                end
            end
            
        end
        
        function self = makeBasisTransformations(self)
            %MAKEBASISTRANSFORMATIONS Creates basis transformation matrices
            
            % Transformation matrix from fully uncoupled to spin-coupled basis
            self.bt21 = zeros(self.Nchannels);
            for row=1:self.Nchannels
                for col=1:self.Nchannels
                    NSpin1Match = self.bv1(col,4)==self.bv2(row,5);
                    NSpin2Match = self.bv1(col,6)==self.bv2(row,6);
                    SpinSumMatch = (self.bv1(col,3)+self.bv1(col,5))==self.bv2(row,4);
                    mLMatch = self.bv1(col,2)==self.bv2(row,2);
                    LMatch = self.bv1(col,1)==self.bv2(row,1);
                    if NSpin1Match && NSpin2Match && SpinSumMatch && mLMatch && LMatch
                        self.bt21(row,col) = ClebschGordan(self.espin,self.espin,self.bv2(row,3),self.bv1(col,3),self.bv1(col,5),self.bv2(row,4));
                    end
                end
            end
            
            % Transformation matrix from fully uncoupled to spin-nucleus coupled basis
            self.bt31 = zeros(self.Nchannels);
            for row=1:self.Nchannels
                for col=1:self.Nchannels
                    SpinSumMatch1 = (self.bv1(col,3)+self.bv1(col,4))==self.bv3(row,4);
                    SpinSumMatch2 = (self.bv1(col,5)+self.bv1(col,6))==self.bv3(row,6);
                    mLMatch = self.bv1(col,2)==self.bv3(row,2);
                    LMatch = self.bv1(col,1)==self.bv3(row,1);
                    if SpinSumMatch1 && SpinSumMatch2 && mLMatch && LMatch
                        tmp1 = ClebschGordan(self.espin,self.nspin(1),self.bv3(row,3),self.bv1(col,3),self.bv1(col,4),self.bv3(row,4));
                        tmp2 = ClebschGordan(self.espin,self.nspin(2),self.bv3(row,5),self.bv1(col,5),self.bv1(col,6),self.bv3(row,6));
                        self.bt31(row,col) = tmp1.*tmp2;
                    end
                end
            end
            
            % Transformation matrix from spin-coupled basis to spin-nucleus coupled basis
            self.bt32 = self.bt31*(self.bt21');
        end
        
        function self = restrict(self,label,dipole)
            %RESTRICT Restricts the basis vectors to an appropriate
            %subspace appropriate to the problem
            %
            %   self = restrict(label,dipole) restricts the basis
            %   vector based on the initial state and whether or not
            %   to include the dipole-dipole interaction.
            %
            %   initLabel is a 3-element vector [L mL mF] where mF is the
            %   total spin projection and dipole is a true/false value.
            %
            %   For dipole = false, subspace restriction keeps states with
            %   the same L, mL, and mF.
            %
            %   For dipole = true, subspace restriction keeps states with
            %   the same mL+mF and L values that have the same parity
            if dipole
                %Restriction to same mL+mF = mJ and L= ' L-2,L,L+2
                matchLabel = [label(1) sum(label(2:3))];
                [~,self.bvint] = self.match(matchLabel,self.bvint,dipole,2:4);
                
                [match1,self.bv1] = self.match(matchLabel,self.bv1,dipole,2:6);
                [match2,self.bv2] = self.match(matchLabel,self.bv2,dipole,[2,4:6]);
                [match3,self.bv3] = self.match(matchLabel,self.bv3,dipole,[2,4,6]);
            else
                %Restriction to same mF AND same L AND mL
                matchLabel = label;
                [~,self.bvint] = self.match(matchLabel,self.bvint,dipole,1,2,3:4);
                
                [match1,self.bv1] = self.match(matchLabel,self.bv1,dipole,1,2,3:6);
                [match2,self.bv2] = self.match(matchLabel,self.bv2,dipole,1,2,4:6);
                [match3,self.bv3] = self.match(matchLabel,self.bv3,dipole,1,2,[4,6]);
            end
            self.bt21 = self.bt21(match2,match1);
            self.bt31 = self.bt31(match3,match1);
            self.bt32 = self.bt31*self.bt21';
            self.L = diag(self.bv1(:,1));
            self.Nchannels = size(self.L,1);
        end
        

        function ops = makeOperators(self,B,initLabel,dipole,refType)
            %MAKEOPERATORS Creates an ALKALIOPERATORS object based on the
            %current basis vectors
            %
            %   ops = makeOperators(B,initLabel,dipole) creates
            %   ALKALIOPERATORS object ops given a magnetic field B (in
            %   Tesla), an initial state label [L mL int1 int2], and
            %   whether or not to use the dipole-dipole interaction dipole
            %
            %   ops = makeOperators(B,initLabel,dipole,'bound') creates
            %   ALKALIOPERATORS object ops given magnetic field B (in
            %   Tesla), an initial state label [L mL mF], and
            %   whether or not to use the dipole-dipole interaction dipole
            Hint0_1 = 0.5*(self.bv3(:,3).*(self.bv3(:,3)+1)-self.espin*(self.espin+1)-self.nspin(1)*(self.nspin(1)+1));
            Hint0_2 = 0.5*(self.bv3(:,5).*(self.bv3(:,5)+1)-self.espin*(self.espin+1)-self.nspin(2)*(self.nspin(2)+1));
            Hint0 = scattconst.h*1e6/scattconst.kb*scattconst.K2A(self.mass)*diag(self.Ahfs(1)*Hint0_1+self.Ahfs(2)*Hint0_2);  %In BV3, as inverse wavenumbers [Angstroms^{-2}]
            Hint0 = self.bt31'*Hint0*self.bt31; %In BV1
            SpinOp = 0.5*diag(self.bv2(:,3).*(self.bv2(:,3)+1)-2*self.espin*(self.espin+1));    %In BV2, no units
            if dipole
                Hdd = (-sqrt(6)*scattconst.alpha.^2.*scattconst.EHartree).*(scattconst.aBohr*1e10)^3*scattconst.K2A(self.mass)/scattconst.kb*DipoleDipole(self.bv1); %In BV1, as inverse wavenumbers*length^3, or Angstroms
            else
                Hdd = zeros(self.Nchannels);
            end
            
            HZ = scattconst.muB*B/scattconst.kb*scattconst.K2A(self.mass)*diag(self.atoms{1}.gS*self.bv1(:,3)+self.atoms{2}.gS*self.bv1(:,5)+...
                 self.atoms{1}.gI*self.bv1(:,4)+self.atoms{2}.gI*self.bv1(:,6));    %In BV1, as inverse wavenumbers [Angstroms^{-2}]
            
            Hint = Hint0+HZ;
            %Determine transformation from basis of internal states in
            %presence of B to fully uncoupled basis
            [~,V1,~,mag1] = self.atoms{1}.hyperfineSolve(B);
            [~,V2,~,mag2] = self.atoms{2}.hyperfineSolve(B);
            BT1int = zeros(self.Nchannels);
            for row=1:self.Nchannels
                for col=1:self.Nchannels
                    if all(self.bv1(row,1:2)==self.bvint(col,1:2))
                        tmp1 = repmat(self.bv1(row,:),size(V1,1),1);
                        tmp2 = repmat(self.bv1(row,:),size(V2,1),1);
                        idx1 = all(tmp1(:,3:4)==mag1,2);
                        idx2 = all(tmp2(:,5:6)==mag2,2);
                        BT1int(row,col) = V1(idx1,self.bvint(col,5))*V2(idx2,self.bvint(col,6));
                    end
                end
            end
            
            BT2int = self.bt21*BT1int;
            %% Change reference level for internal Hamiltonian
            Hint = BT1int'*Hint*BT1int;
            if nargin<5 || strcmpi(refType,'scatt')
                initIdx = self.findstate(self.bvint(:,[1,2,5,6]),initLabel);
            elseif strcmpi(refType,'bound')
                [~,initIdx] = min(diag(Hint));
            else
                error('Energy reference type not recognized');
            end
            refE = Hint(initIdx,initIdx);
            
            %% Change operators to BV2 (total electron spin)
            Hint0 = self.bt21*Hint0*self.bt21';
            Hdd = self.bt21*Hdd*self.bt21';
            Hint = BT2int*Hint*BT2int';
            
            ops = alkalioperators(self.pair,self.scale,self.L,SpinOp,Hdd,Hint0,Hint,refE,BT2int);
            
            
        end
        
        
    end
    
    methods(Static)
        function idx = findstate(bv,vec)
            %FINDSTATE Returns the index corresponding to a vector in a
            %set of basis vectors
            %
            %   idx = atompairbasis.findstate(bv,vec) finds the index idx
            %   corresponding to vec in bv
            idx = all(bv==repmat(vec,size(bv,1),1),2);
        end
        
        function [match,bvout] = match(matchLabel,bv,dipole,varargin)
            % MATCH Finds states in BV that have the same values using
            % matchLabel
            %   [match,bvout]=match(matchLabel,bv,dipole,i1,i2,...,iN)
            %   finds all states in bv with the same quantum numbers
            %   defined by matchLabel and whether or not to include the
            %   dipole-dipole interaction dipole.  The indices in bv to
            %   look at are defined by i1,i2,...,iN.  If any of iN are
            %   vectors with more than one element, the matched quantity is
            %   summed.  So, for instance, matchLabel=[0,0,5] and i1=1,
            %   i2=2, and i3=3:4 looks for rows in bv where the first
            %   element is equal to 0, the second is equal to 0, and the
            %   sum of the third and fourth elements is equal to 5.  Vector
            %   M is a logical vector that selects the matching elements,
            %   and BVout is the restricted version of BV that satisfies
            %   the matching condition.
            %
            %   [match,bvout]=match(matchLabel,bv,dipole,idx) matches
            %   quantum numbers according to dipole-dipole selection rules.
            %   The first element in matchLabel is assumed to be the value
            %   of L, and bv rows with L values that have the same parity.
            %   idx is the vector of indices defining all the angular
            %   momentum projections, such that mJ=sum(bv(:,idx),2)
            if dipole
                args = varargin;
                L = matchLabel(1);
                matchLabel = matchLabel(2:end);
            else
                args = varargin;
            end
            
            if numel(args)~=numel(matchLabel)
                error('Matching labels have different lengths');
            end
            
            match = false(size(bv,1),numel(args));
            for nn=1:numel(args)
                match(:,nn) = repmat(matchLabel(nn),size(bv,1),1)==sum(bv(:,args{nn}),2);
            end
            
            match = all(match,2);
            
            if dipole
                match = match & (mod(bv(:,1)-L,2)==0);
            end
            
            bvout = bv(match,:);
        end
    end
    
    
end