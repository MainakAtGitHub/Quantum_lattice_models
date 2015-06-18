function [Impparameters, ImpVector] = process_impurity_tb_tom( filename, kz,swap)
fid = fopen(filename);
if nargin < 2
    kz=NaN
end;
if nargin <3
    swap= false;
end;
% read string in first line
oldformat=false;
if oldformat
    %[s]= fscanf(fid, '%s', 10)
    [number]= fscanf(fid, ' %d', 1);
else
    s='';

    while ~(strcmp('substitutions:',s))
        [s]= fscanf(fid, '%s', 1)
    end;
    % read number of matrices and number of orbitals, ignoring the first
    % lines
[number]= fscanf(fid, ' %d ', 1)
end;
[nOrb]= fscanf(fid, '%d', 1)
readstring1= ['\t(%g,%g)\t'];
%for o=1:nOrb
 %   readstring1=[readstring1,' (%g,%g)'];
%end;
number1=0
for n=1:number
    rrp= fscanf(fid, '\t%d', 6);
    line1=fscanf(fid,readstring1,2*nOrb^2);
    line2=reshape(line1,2,nOrb^2)
    line3=line2(1,:)+1i*line2(2,:);
    matrix=reshape(line3,nOrb,nOrb);
   % lv=[rrp(1)-rrp(4),rrp(2)-rrp(5),rrp(3)-rrp(6)];
    lv=[rrp(1:6)']:
    if ~isnan(kz)
    if (lv(3)==0) && (lv(6)==0) 
        number1=number1+1;
        %latticeVector(number1,:)=[rrp(1)-rrp(4),rrp(2)-rrp(5),rrp(3)-rrp(6)];
                ImpVector(number1,:)=lv;
                % do a transpose due to convention differences!
        Impparameters(:,:,number1)=matrix;
    else
       number1=number1+1;
       lv(3)=0;
       lv(6)=0;
       ImpVector(number1,:)=lv;
       Impparameters(:,:,number1)=matrix*exp(1i*(rrp(3)-rrp(6))*kz*pi);
    end;
    end;
end;
    % only consider real part?
im=max(abs(imag(Impparameters(:))));
if im < 1e-8
    disp(['Only taking real part, imaginary part is small:',num2str(im)]);
    Impparameters=real(Impparameters);
end;
% swap Fe(1) <-> Fe(2)
if swap
    disp('Swapping Fe(1) and Fe(2)');
    Impparameters=swap_matrix(Impparameters);
end;
save([filename,'_conv_z',num2str(kz),'.mat'],'Impparameters','ImpVector');
% number2=0;
% for orb1=1:nOrb
%     for orb2=1:nOrb
%         for r=1:number1
%             if abs(Impparameters(orb1,orb2,r))>0
%                 number2=number2+1;
%                 tb(number2,:)=[ImpVector(r,:),orb1,orb2,Impparameters(orb1,orb2,r)];
%             end;
%         end
%     end;
% end
% dlmwrite([filename,'_conv_z',num2str(kz),'.csv'],tb,'precision',10)
% end

