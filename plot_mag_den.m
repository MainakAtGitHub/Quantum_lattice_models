
% imagesc(reshape(nUp(1:1:1*12^2)-nDown(1:1:1*12^2),[12,12]));
% surf(reshape(nUp(1:1:1*16^2)-nDown(1:1:1*16^2),[16,16]));
% surf(reshape(nUp(1:1:1*16^2)+nDown(1:1:1*16^2),[16,16]));
% 
% imagesc(reshape(nUp(1:2:2*11*11)+nUp(2:2:2*11*11)-nDown(1:2:2*11*11)-nDown(2:2:2*11*11),[11,11]))
% 
% plot([1, 2, 3, 4,1.25,1.5,1.75,2.25],[-.3461, .0789, .4171, .5257,-0.239884425593240,-0.133634425637637,-0.027384425637727,0.183027619935147],'*')
% 
% surf(reshape(nUp(1:1:nOrbitals*N^2)-nDown(1:1:nOrbitals*N^2),[N,N]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Above is some test code,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ignore
function f = plot_mag_den(inputfile,orbNo,dislocation_length)

if nargin<3
    dislocation_length=false;
end

read_input_file=inputfile;
read_input;
read_input_file

load(BdGfileName,'-mat');
if dislocation_length>0
    nUpAux=zeros(N^2,1);
    nDownAux=zeros(N^2,1);
    itr=0;
    for i=1:N
        for j=1:N
            if ~and(j==(N+1)/2,any((N+1)/2-(dislocation_length-1)/2:(N+1)/2+(dislocation_length-1)/2==i))
                itr=itr+1;
                nUpAux(N*(i-1)+j)=nUp(itr);
                nDownAux(N*(i-1)+j)=nDown(itr);
            else
                nUpAux(N*(i-1)+j)=0;
                nDownAux(N*(i-1)+j)=0;
            end
        end
    end
    nUp=nUpAux;
    nDown=nDownAux;
end

nUpAllOrb = zeros(length(nUp)/nOrbitals,1);
nDownAllOrb = zeros(length(nDown)/nOrbitals,1);

if nargin < 2
    orbNo=0;
end
if orbNo==0
    for orbNo = 1:nOrbitals
        nUpAllOrb = nUpAllOrb + nUp(orbNo:nOrbitals:nOrbitals*N^2);
        nDownAllOrb = nDownAllOrb + nDown(orbNo:nOrbitals:nOrbitals*N^2);
    end 
else
        nUpAllOrb = nUpAllOrb + nUp(orbNo:nOrbitals:nOrbitals*N^2);
        nDownAllOrb = nDownAllOrb + nDown(orbNo:nOrbitals:nOrbitals*N^2);
    
end

mag = nUpAllOrb-nDownAllOrb;

stag_mag = zeros(N);
temp_reshaped_mag = reshape(mag,[N,N])';
for ii=1:N
    for jj=1:N
        stag_mag(ii,jj)=(-1)^(ii+jj)*temp_reshaped_mag(ii,jj);
    end
end

den = nUpAllOrb+nDownAllOrb;

[filepath,name,ext]=fileparts(BdGfileName);
figure;
% surf(reshape(mag,[N,N])');
imagesc(reshape(mag,[N,N]));
title('Magnetization')
axis square;
colorbar;

print_pdf([filepath,filesep,name,'_magnetization.pdf']);
% saveas(gcf, 'Magnetization');
% saveas(gcf, 'Magnetization.jpg');

figure;
% surf(stag_mag);
imagesc(transpose(stag_mag));
% imagesc(reshape(stag_mag,[N,N])');
title('Staggered Magnetization')
axis square;
colorbar;
print_pdf([filepath,filesep,name,'_staggered_magnetization.pdf']);
% saveas(gcf, 'Staggered Magnetization');
% saveas(gcf, 'Staggered_Magnetization.jpg');

figure;
% surf(reshape(den,[N,N])');
imagesc(reshape(den,[N,N]));
title('Density')
axis square;
colorbar;
print_pdf([filepath,filesep,name,'_density.pdf']);
% saveas(gcf, 'Density');
% saveas(gcf, 'Density.jpg');