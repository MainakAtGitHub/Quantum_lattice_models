function r=write_complex_matrix(flnm,data)
if ~isempty(flnm)
    fid = fopen(flnm,'wt');
    % make sure the file is not empty
    finfo = dir(flnm);
    fsize = finfo.bytes;
else
    r=[];
end;
fprintf('Writing %u values to file %s.\n', size(data,2), flnm)
for n=1:size(data,2)
    for m=1:size(data,1)
        if ~isempty(flnm)
            fprintf(fid,'(%e,%e)',real(data(m,n)),imag(data(m,n)));
        else
            r=[r,sprintf('(%e,%e)',real(data(m,n)),imag(data(m,n)))];
        end;
    end
    if ~isempty(flnm)
        fprintf(fid,'\n');
    else
        r=[r,sprintf('\n')];
    end;
end
if ~isempty(flnm)
    fclose(fid);
end