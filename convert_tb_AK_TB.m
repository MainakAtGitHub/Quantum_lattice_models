function [TBparameters,latticeVector]=convert_tb_AK_TB(tbfile)
tb=load(tbfile,'-ascii');
%cut=0;
nOrb=max(max(tb(:,4:5)))
sztb=size(tb);
TBparameters=[];
latticeVector=[];
for n=1:sztb(1)
    if ~isempty(TBparameters)
       % n
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
    TBparameters(tb(n,4),tb(n,5),rowColIdx)=TBparameters(tb(n,4),tb(n,5),rowColIdx)+tb(n,6);
end;
sublattice=input('enter sublattice 0,1,-1: ');
% check for CC
%tb1=permute(TBparameters,[2,1,3]);
% if sum(abs(tb1(:)-TBparameters(:)))>1e-3
%     disp('adding cc!');
%     tbtmp=TBparameters;
%     for i=1:size(TBparameters,3)
%         % figure out r and -r
%         lattvec=-repmat(latticeVector(i,:),size(latticeVector,1),1);
%         negindex=find(sum(latticeVector==lattvec,2)==size(latticeVector,2))        
%         TBparameters(:,:,i)=0.5*(tbtmp(:,:,i)+tb1(:,:,negindex));%-diag(diag(tbtmp(:,:,i)));
%     end
% end
% write in the format of Tom Berlijn
[latticeVector,ind]=sortrows(latticeVector);
TBparameters=TBparameters(:,:,ind);
outfile=[tbfile,'conv.out']
    fid = fopen(outfile,'wt');
    % make sure the file is not empty
    finfo = dir(outfile);
    fsize = finfo.bytes;
    % print a newline
    fprintf(fid,'\n');
    % print the geometry
    r=sprintf('%i %i\n',size(TBparameters,3),size(TBparameters,2));
    fprintf(fid,'%s',r);
    for i=1:size(TBparameters,3)
        r=sprintf('%i %i %i\n',latticeVector(i,:));
        for orb1=1:size(TBparameters,2)
           for orb2=1:size(TBparameters,1) 
               r=[r,sprintf('(%e,%e)',real(TBparameters(orb1,orb2,i)),imag(TBparameters(orb1,orb2,i)))];
           end
           r=[r,sprintf('\n')];
        end
        fprintf(fid,'%s',r);
    end;
    fclose(fid);

    

%save([tbfile,'.mat'],'TBparameters','latticeVector','sublattice');
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
