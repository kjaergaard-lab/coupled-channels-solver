classdef ScatteringMatrix < matlab.mixin.Copyable
    %ScatteringMatrix Provides a convenient interface for common tasks
    %involving scattering matrix calculations
    %
    %   The main convenience is how different elements in the S or T
    %   matrices are accessed.  See help for ScatteringMatrix.subsref for
    %   details
    properties
        mass        %The reduced mass of the system in question
        E           %The energies (in uK)
        B           %The magnetic fields (in G)
        k           %The calculated wavenumbers
        
        bv          %The basis vectors (internal)
        symmetry    %The symmetry +/-1 for bosons/fermions or 0 for others
        Sfull       %The full S-matrix (Nch x Nch x Nruns)
        Tfull       %The full T-matrix (Nch x Nch x Nruns)
        targetIndex %The index of bv indicating the entrance channel     
    end
    
    methods
        function self = ScatteringMatrix(Sin,BVin,symIn,targetIndex)
            %SCATTERINGMATRIX Constructs an instance of the
            %ScatteringMatrix class
            %
            %   self = ScatteringMatrix(Sin,BVin,symIn,targetIndex)
            %   constructs an instance of ScatteringMatrix with the input
            %   S-matrix Sin, basis vectors BVin, symmetry symIn=0,-1,1,
            %   and entrance channel specified by targetIndex corresponding
            %   to the entrance channel row in BVin
            self.bv = BVin;
            self.symmetry = symIn;
            self.targetIndex = targetIndex;
            self.Sfull = Sin;
            self.Tfull = Sin-repmat(eye(size(Sin(:,:,1))),[1,1,size(Sin,3)]); 
        end
        
        function set.E(self,E)
            %set.E Sets the energy (in uK) and calculates the wavenumber
            tmp = E(:)*scattconst.kb/1e6;
            self.E = E(:);
            self.k = sqrt(2*self.mass/scattconst.hbar^2*tmp(:));  %#ok
        end
        
        function cs = crossSec(self,varargin)
            %CROSSSEC Calculates and returns various cross sections.
            %
            %   cs = crossSec() calculates the elastic cross section
            %   starting and ending in the entrance channel
            %
            %   cs = crossSec(finalState) where finalState is a 1x2 vector
            %   calculates the cross section for collisions starting in the
            %   entrance channel states and ending in the final internal
            %   states given by finalState
            %
            %   cs = crossSec('total') calculates the total cross section
            %   for atoms starting in the entrance channel
            %
            %   cs = crossSec('inelastic') calculates the total inelastic
            %   cross section (total - elastic)
            %
            %   cs = crossSec(initState,finalState) calculates the cross
            %   section for atoms starting in initState and ending in
            %   finalState.  Both initState and finalState must be 1x2
            %   vectors specifying states in the internal basis
            if nargin == 1
                initState = self.bv(self.targetIndex,5:6);
                finalState = initState;
            elseif nargin == 2
                if ~ischar(varargin{1})
                    initState = self.bv(self.targetIndex,5:6);
                    finalState = varargin{1};
                else
                    [~,ia] = unique(self.bv(:,5:6),'stable','rows');
                    intState = self.bv(ia,:);
                    cs = zeros(numel(self.E),1);
                    for nn=1:size(intState,1)
                        cs = cs+self.crossSec(intState(nn,5:6));
                    end
                    
                    if strcmpi(varargin{1},'total')
                        return
                    elseif strcmpi(varargin{1},'inelastic')
                        cs = cs-self.crossSec();
                        return;
                    end
                end
            elseif nargin == 3
                initState = varargin{1};
                finalState = varargin{2};
            end
            cs = self.partialCrossSec(self.Tfull,self.bv(:,[1,2,5,6]),self.k,initState,finalState,self.symmetry);
        end
        
        function self = plotCrossSections(self)
            %PLOTCROSSSECTIONS Plots all partial cross sections starting
            %from the entrance state
            [~,ia] = unique(self.bv(:,5:6),'stable','rows');
            intState = self.bv(ia,:);
            jj = 1;
            clf;
            for nn=1:size(intState,1)
                cs = self.crossSec(intState(nn,5:6));
                if cs(end)~=0
                    if numel(self.E)>1
                        plot(self.E,cs,'.-');
                    else
                        plot(self.B,cs,'.-');
                    end
                    hold on;
                    str{jj} = sprintf('%d-%d',intState(nn,5:6)); %#ok<AGROW>
                    jj=jj+1;
                end
            end
            if numel(self.E) > 1
                plot(self.E,self.crossSec('inelastic'),'k--','LineWidth',1.5);
            else
                plot(self.B,self.crossSec('inelastic'),'k--','LineWidth',1.5);
            end
            hold off;
            legend(str);
            fs = 10;
            xlabel('Collision Energy [uK]','FontSize',fs);
            ylabel('Cross section [m^2]','FontSize',fs);
        end

        function R = collisionRate(self,varargin)
            %COLLISIONRATE Calculates and returns various collision rate
            %coefficients
            %
            %   R = collisionRate() calculates the elastic collision rate
            %   coefficient starting and ending in the entrance channel
            %
            %   cs = collisionRate(finalState) where finalState is a 1x2
            %   vector calculates the collision rate coefficient for
            %   collisions starting in the entrance channel states and
            %   ending in the final internal states given by finalState
            %
            %   cs = collisionRate('total') calculates the total collision
            %   rate coefficient for atoms starting in the entrance channel
            %
            %   cs = collisionRate('inelastic') calculates the total
            %   inelastic collision rate coefficient (total - elastic)
            %
            %   cs = collisionRate(initState,finalState) calculates the
            %   collision rate coefficient for atoms starting in initState
            %   and ending in finalState.  Both initState and finalState
            %   must be 1x2 vectors specifying states in the internal basis

            cs = self.crossSec(varargin{:});
            R = cs.*scattconst.hbar.*self.k/self.mass;
        end

        function self = plotCollisionRates(self)
            %PLOTCOLLISIONRATES Plots all partial collision rate
            %coefficients starting from the entrance state
            [~,ia] = unique(self.bv(:,5:6),'stable','rows');
            intState = self.bv(ia,:);
            jj = 1;
            for nn=1:size(intState,1)
                R = self.collisionRate(intState(nn,5:6));
                if R(end)~=0
                    if numel(self.E)>1
                        plot(self.E,R,'.-');
                    else
                        plot(self.B,R,'.-');
                    end
                    hold on;
                    str{jj} = sprintf('%d-%d',intState(nn,5:6)); %#ok<AGROW>
                    jj=jj+1;
                end
            end
            if numel(self.E) > 1
                plot(self.E,self.collisionRate('inelastic'),'k--','LineWidth',1.5);
            else
                plot(self.B,self.collisionRate('inelastic'),'k--','LineWidth',1.5);
            end
            hold off;
            legend(str);
            fs = 10;
            xlabel('Collision Energy [uK]','FontSize',fs);
            ylabel('Two-body rate coefficient [m^3/s]','FontSize',fs);
        end
        
        function B = subsref(self,S)
            %SUBSREF Controls how subscript referencing works for
            %ScatteringMatrix class
            %
            %   Overloading this function allows for a convenient method by
            %   which different elements in the S or T matrices can be
            %   accessed.
            %
            %   General usage for an instance S of ScatteringMatrix is
            %   S(subscripts).(Property).  If (Property) is something other
            %   than S or T, the default behaviour applies, so one can
            %   access bv or Sfull or Tfull in a sensible manner.
            %
            %   If one uses S.S or S.T, it returns the elastic element of the
            %   S or T matrix, which is entrance channel to entrance
            %   channel.
            %
            %   If one uses S(idx).S or S(idx).T then the element that is
            %   accessed is actually S(entrance channel index,idx)
            %
            %   If one uses S(L,mL,int1,int2).S then the element that is
            %   accessed is entrance channel to the channel specified by
            %   the 4-element vector [L,mL,int1,int2].
            %
            %   If one uses S(idx1,idx2).S then the element that is
            %   accessed is the one going from idx1 to idx2 where these are
            %   indices in the basis vector bv
            %
            %   If one uses S(L1,mL1,int1_1,int2_1,L2,mL2,int1_2,int2_2).S
            %   then the element that is accessed is the one starting from
            %   the 4-element vector [L1,mL1,int1_1,int2_1] and ending with
            %   [L2,mL2,int1_2,int2_2]
            if numel(self) > 1
                B = builtin('subsref',self,S);
            end
            idx1 = self.targetIndex;
            idx2 = idx1;
            for nn=1:numel(S)
                switch S(nn).type
                    case '.'
                        switch S(nn).subs
                            case 'S'
                                B = squeeze(self.Sfull(idx1,idx2,:));
                            case 'T'
                                B = squeeze(self.Tfull(idx1,idx2,:));
                            otherwise
                                B = builtin('subsref',self,S);
                                return;
                        end
                        
                    case '()'
                        if nn==1
                            if numel(S(nn).subs)==1
                                idx2 = S(nn).subs{1};
                            elseif numel(S(nn).subs)==2
                                idx1 = S(nn).subs{1};
                                idx2 = S(nn).subs{2};
                            elseif numel(S(nn).subs)==4
                                idx2 = find(all(self.bv(:,[1,2,5,6])==repmat(cell2mat(S(nn).subs),size(self.bv,1),1),2));
                            elseif numel(S(nn).subs)==8
                                v = cell2mat(S(nn).subs);
                                v1 = v(1:4);v2 = v(5:8);
                                idx1 = find(all(self.bv(:,[1,2,5,6])==repmat(v1,size(self.bv,1),1),2));
                                idx2 = find(all(self.bv(:,[1,2,5,6])==repmat(v2,size(self.bv,1),1),2));
                            else
                                error('Index type not supported');
                            end
                        else
                            tmp = B;
                            B = subsref(tmp,S(nn));
                        end
                        
                    otherwise
                        error('Index type not supported');
                end
            end
        end
    end
    
    methods(Static)
       function cs = partialCrossSec(T,bv,k,initState,finalState,identical)
           % partialCrossSec Calculates partial cross sections
           %   cross_sec=partialCrossSec(T,bv,k,initState,finalState) with T the
           %   full T-matrix, bv the basis vectors of T, k the wavenumber, and
           %   initState and finalState are denoted by the [Int1,Int2] state labels
           
           InitIdx = atompairbasis.findstate(bv(:,3:4),initState);
           FinalIdx = atompairbasis.findstate(bv(:,3:4),finalState);
           
           InLRange=unique(bv(InitIdx,1));
           InMLRange=unique(bv(InitIdx,2));
           OutLRange=unique(bv(FinalIdx,1));
           OutMLRange=unique(bv(FinalIdx,2));
           
           cs=zeros(size(T,3),1);
           for n1=1:numel(InLRange)
               for n2=1:numel(InLRange)
                   for n3=1:numel(OutLRange)
                       for n4=1:numel(OutMLRange)
                           for n5=1:numel(InMLRange)
                               for n6=1:numel(InMLRange)
                                   L1=InLRange(n1);
                                   L2=InLRange(n2);
                                   L=OutLRange(n3);
                                   mL=OutMLRange(n4);
                                   mL1=InMLRange(n5);
                                   mL2=InMLRange(n6);
                                   idx1=atompairbasis.findstate(bv,[L1,mL1,initState]);
                                   idx2=atompairbasis.findstate(bv,[L2,mL2,initState]);
                                   idx3=atompairbasis.findstate(bv,[L,mL,finalState]);
                                   cs=cs+squeeze(1i^(L2-L1).*sqrt((2*L1+1).*(2*L2+1)).*conj(T(idx3,idx1,:)).*T(idx3,idx2,:));
                               end
                           end
                       end
                   end
               end
           end
           if nargin <6
               identical = 0;
           end
           
           cs = (1+identical*(finalState(1)==finalState(2)))*pi./k.^2.*real(cs);
       end
    end
    
    
end