function [TBparameters, latticeVector, sublattice] = process_tb_tom( filename, kz)
fid = fopen(filename);
if nargin < 2
    kz=NaN
end;
% read string in first line
oldformat=false;
if oldformat
    [s]= fscanf(fid, '%s\t', 13)
    [number]= fscanf(fid, ' %d', 1);
else
    % read number of matrices and number of orbitals
[number]= fscanf(fid, ' %d', 1);
end;
[nOrb]= fscanf(fid, '%d', 1);
if oldformat
[s]= fscanf(fid, '%s', 1)
[s]= fscanf(fid, '%s', 1)
[s]= fscanf(fid, '%s', 1)
[s]= fscanf(fid, '%s', 1)
[s]= fscanf(fid, '%s', 1)
end;
sublattice=-1;
readstring1= ['\t(%g,%g)\t'];
%for o=1:nOrb
 %   readstring1=[readstring1,' (%g,%g)'];
%end;
number1=0
for n=1:number
    if oldformat
        rrp= fscanf(fid, '\t%d', 6);
    else
             rrp= fscanf(fid, '\t%d', 3);   
    end
    line1=fscanf(fid,readstring1,2*nOrb^2);
    line2=reshape(line1,2,nOrb^2);
    line3=line2(1,:)+1i*line2(2,:);
    matrix=reshape(line3,nOrb,nOrb);
   % lv=[rrp(1)-rrp(4),rrp(2)-rrp(5),rrp(3)-rrp(6)];
    lv=[rrp(1),rrp(2),rrp(3)];
    if isnan(kz)
    if lv(3)==0
        number1=number1+1;
        %latticeVector(number1,:)=[rrp(1)-rrp(4),rrp(2)-rrp(5),rrp(3)-rrp(6)];
                latticeVector(number1,:)=[rrp(1),rrp(2),rrp(3)];
                % do a transpose due to convention differences!
        %TBparameters(:,:,number1)=matrix';
                TBparameters(:,:,number1)=matrix;
    end;
    else
       number1=number1+1;
       latticeVector(number1,:)=[rrp(1),rrp(2),0];
%       TBparameters(:,:,number1)=matrix'*exp(1i*rrp(3)*kz*pi);
              TBparameters(:,:,number1)=matrix*exp(1i*rrp(3)*kz*pi);
    end;
end;
    % only consider real part?
im=max(abs(imag(TBparameters(:))));
if im < 1e-8
    TBparameters=real(TBparameters);
end;
save([filename,'_conv_z',num2str(kz),'a.mat'],'TBparameters','latticeVector','sublattice');
number2=0;
for orb1=1:nOrb
    for orb2=1:nOrb
        for r=1:number1
            if abs(TBparameters(orb1,orb2,r))>0
                number2=number2+1;
                tb(number2,:)=[-latticeVector(r,:),orb1,orb2,TBparameters(orb1,orb2,r)];
            end;
        end
    end;
end
dlmwrite([filename,'_conv_z',num2str(kz),'a.csv'],tb,'precision',10)
end

