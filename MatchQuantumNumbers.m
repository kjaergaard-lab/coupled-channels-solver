function [M,BVout]=MatchQuantumNumbers(MatchLabel,BV,varargin)
% MatchQuantumNumbers Finds states in BV that have the same values using MatchLabel
%   [M,BVout]=MatchQuantumNumbers(MatchLabel,BV,i1,i2,...,iN) finds all states
%   in BV with the same quantum numbers defined by MatchLabel.  The indices in 
%   BV to look at are defined by i1,i2,...,iN.  If any of iN are vectors
%   with more than one element, the matched quantity is summed.  So, for
%   instance, MatchLabel=[0,0,5] and i1=1, i2=2, and i3=3:4 looks for rows
%   in BV where the first element is equal to 0, the second is equal to 0,
%   and the sum of the third and fourth elements is equal to 5.  Vector M
%   is a logical vector that selects the matching elements, and BVout is
%   the restricted version of BV that satisfies the matching condition.
%
%   [M,BVout]=MatchQuantumNumbers(MatchLabel,BV,idx,'Dipole') matches
%   quantum numbers according to dipole-dipole selection rules.  The first
%   element in MatchLabel is assumed to be the value of L, and BV rows with
%   L'=L-2,L,L+2 are kept.  idx is the vector of indices defining all the
%   angular momentum projections, such that mJ=sum(BV(:,idx),2)

if ischar(varargin{end}) && strcmpi(varargin{end},'Dipole'),
    DipoleFlag=1;
    ArgList=varargin(1:end-1);
    LLabel=MatchLabel(1);
    MatchLabel=MatchLabel(2:end);
else
    DipoleFlag=0;
    ArgList=varargin;
end;

if numel(ArgList)~=numel(MatchLabel),
    error('Matching labels have different lengths!');
end;

M=false(size(BV,1),numel(ArgList));
for nn=1:numel(ArgList),
    idx=ArgList{nn};
    M(:,nn)=repmat(MatchLabel(nn),size(BV,1),1)==sum(BV(:,idx),2);
end;

M=all(M,2);
if DipoleFlag,
    M2=(BV(:,1)==LLabel) | (BV(:,1)==(LLabel-2)) | (BV(:,1)==(LLabel+2));
    M=M & M2;
end;
    
BVout=BV(M,:);
