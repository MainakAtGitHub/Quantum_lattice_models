function c=convert_delta_AK(filename,N,nOrb)
if nargin <2
    N=9;
end
if nargin <3
    nOrb=10;
end;
load(filename);
center=ceil(N/2);
filehandle=fopen([filename,'conv'],'w');
for i=1:11
    for j=1:11
        nx=center-6+i;
        ny=center-6+j;
        if (nx>0) && (ny >0) && (nx < N+1) && (ny < N+1)
            [nx ny]
            [iRange, jRange]=find_lattice_translation_index(N,nOrb,[ center center],[nx ny]);
            spatialgap=delta(iRange,jRange);
        else
            spatialgap=zeros(nOrb,nOrb);
        end;
        strng=[];
        for index=1:numel(spatialgap)
            strng=[strng,'(',num2str(real(spatialgap(index))),',',num2str(imag(spatialgap(index))),')'];
        end;
        fprintf(filehandle,'%s\n',strng);
    end;
end;
fprintf(filehandle,'%s\n',[]);
fclose(filehandle);