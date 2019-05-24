function [TBparameters,latticeVector]=convert_tb_AK_PC(tbfile)
try
tb=load(tbfile,'-ascii');
catch
    load(tbfile,'tb','-mat');
end;
%cut=0;
nOrb=max(max(tb(:,4:5)))
sztb=size(tb);
TBparameters=[];
latticeVector=[];
for n=1:sztb(1)
    if ~isempty(TBparameters)
        n
        [ rowColIdx ] = findRowOrColumnInMatrix( latticeVector,tb(n,1:3));
    else
        rowColIdx=[];
    end;
    if isempty(rowColIdx)
        % add the lattice vector
        latticeVector=[latticeVector;tb(n,1:3)];
        rowColIdx=size(latticeVector,1);   
        TBparameters(nOrb,nOrb,rowColIdx)=0;
    end;
    try
        TBparameters(tb(n,4),tb(n,5),rowColIdx)=TBparameters(tb(n,4),tb(n,5),rowColIdx)+tb(n,6)+1i*tb(n,7);
    catch
        TBparameters(tb(n,4),tb(n,5),rowColIdx)=TBparameters(tb(n,4),tb(n,5),rowColIdx)+tb(n,6);
    end;
end;
sublattice=input('enter sublattice 0,1,-1: ');
% check for CC
tb1=permute(TBparameters,[2,1,3]);
if sum(abs(tb1(:)))-sum(abs(TBparameters(:)))>1e-3
    disp('adding cc!');
    tbtmp=TBparameters;
    for i=1:size(TBparameters,3)
        % figure out r and -r
        lattvec=-repmat(latticeVector(i,:),size(latticeVector,1),1);
        negindex=find(sum(latticeVector==lattvec,2)==size(latticeVector,2))        
        TBparameters(:,:,i)=0.5*(tbtmp(:,:,i)+tb1(:,:,negindex));%-diag(diag(tbtmp(:,:,i)));
    end
end
save([tbfile,'.mat'],'TBparameters','latticeVector','sublattice');
%tb=[];
%sztb=size(TBparameters);
%for mu=1:sztb(1)
%    for nu=1:sztb(2)
%        for r=1:sztb(3)
%            if abs(TBparameters(mu,nu,r))>cut
%                tb=[tb;latticeVector(r,:),0,mu,nu,TBparameters(mu,nu,r)];
%            end;
%        end
%    end
%end
