function [TBparameters,latticeVector]=convert_tb_AK_PC(tbfile,sublattice)
% read in a standart tight binding file (similar to output from FPLO,
% Wannier90 and convert to the matlab file as introtuced by Peayush
% second argument sets the "sublattice" value that is required for 10
% orbital calculations on FeSC, if ommitted, the program will ask for user
% input

% also allow to load matlab files containing the matrix
try
    tb=load(tbfile,'-ascii');
catch
    load(tbfile,'tb','-mat');
end
%cut=0;
% read off the number of orbitals
nOrb=max(max(tb(:,4:5)))
% obtain the number of hoppings
sztb=size(tb);
% initialize variables to be empty
TBparameters=[];
latticeVector=[];
% iterate over hopping elements and store them into structure of matrices
% and the corresponding real space vectors
for n=1:sztb(1)
    % check whether lattice vector is already in list
    if ~isempty(TBparameters)
        n
        [ rowColIdx ] = findRowOrColumnInMatrix( latticeVector,tb(n,1:3));
    else
        rowColIdx=[];
    end
    if isempty(rowColIdx)
        % add the lattice vector
        latticeVector=[latticeVector;tb(n,1:3)];
        rowColIdx=size(latticeVector,1);   
        TBparameters(nOrb,nOrb,rowColIdx)=0;
    end
    try
        % case for complex hoppings as input
        TBparameters(tb(n,4),tb(n,5),rowColIdx)=TBparameters(tb(n,4),tb(n,5),rowColIdx)+tb(n,6)+1i*tb(n,7);
    catch
        % otherwise hoppings are real
        TBparameters(tb(n,4),tb(n,5),rowColIdx)=TBparameters(tb(n,4),tb(n,5),rowColIdx)+tb(n,6);
    end
end
if nargin < 2
    % special cases for sublattice calculations (1Fe cell vs. 2Fe cell):
    sublattice=input('enter sublattice 0,1,-1: ');
end
% check for complex conjugate, if not present, add the corresponding terms
% which sometimes are removed from the input tight-binding files
tb1=permute(TBparameters,[2,1,3]);
if sum(abs(tb1(:)))-sum(abs(TBparameters(:)))>1e-3
    % display a message when adding complex conjugate
    disp('adding cc!');
    tbtmp=TBparameters;
    for i=1:size(TBparameters,3)
        % figure out lattice vector r and its negative -r
        lattvec=-repmat(latticeVector(i,:),size(latticeVector,1),1);
        negindex=find(sum(latticeVector==lattvec,2)==size(latticeVector,2))        
        TBparameters(:,:,i)=0.5*(tbtmp(:,:,i)+tb1(:,:,negindex));%-diag(diag(tbtmp(:,:,i)));
    end
end
% save the file into a new file (matlab format)
save([tbfile,'.mat'],'TBparameters','latticeVector','sublattice');