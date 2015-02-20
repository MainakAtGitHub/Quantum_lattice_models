function r=write_complex_matrix(flnm,data)
fid = fopen(flnm,'wt');
% make sure the file is not empty
finfo = dir(flnm);
fsize = finfo.bytes;
fprintf('Writing %u values to file %s.\n', size(data,2), flnm)
for n=1:size(data,2)
    for m=1:size(data,1)
        fprintf(fid,'(%e,%e)',real(data(m,n)),imag(data(m,n)));
    end
    fprintf(fid,'\n');
end
fclose(fid);