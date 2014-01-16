function c=convert_delta_AK(filename,nOrb)
%if nargin <2
%    N=9;
%end
if nargin <2
    nOrb=10;
end;
load(filename);
szdelta=size(delta);
N=sqrt(szdelta(1)/10);
center=ceil(N/2);
filehandle=fopen([filename,'conv'],'w');
N1=N;
for i=1:N1
    for j=1:N1
        nx=i;
        ny=j;
        if (nx>0) && (ny >0) && (nx < N+1) && (ny < N+1)
            [nx ny]
            [iRange, jRange]=find_lattice_translation_index(N,nOrb,[nx ny],[center center]);
            spatialgap=delta([iRange(6:10),iRange(1:5)],[jRange(6:10),jRange(1:5)])';
        else
            spatialgap=zeros(nOrb,nOrb);
        end;
        strng=[];
        for index=1:numel(spatialgap)
%            strng=[strng,'(',num2str(real(spatialgap(index))),',',num2str(imag(spatialgap(index))),')'];
                        strng=[strng,'(',num2str(real(spatialgap(index))),',','0)'];

        end;
        fprintf(filehandle,'%s\n',strng);
    end;
end;
fprintf(filehandle,'%s\n',[]);
fclose(filehandle);
% also write out the corresponding real space grid (for calculatiion
% purpose)
[Nx,Ny]=meshgrid(-(center-1):center-1);
Nz=0*Nx;
gridN=[Nx(:),Ny(:),Nz(:)];
csvwrite(['real_grid_N',num2str(N),'.csv'],gridN);
