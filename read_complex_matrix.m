function r=read_complex_matrix(flnm,rows,columns,skip)
fid = fopen(flnm,'rt');
% make sure the file is not empty
finfo = dir(flnm);
fsize = finfo.bytes;
if ~exist('rows','var')
    rows=1;
end;
if ~exist('skip','var')
    skip=0;
end;
if ~exist('columns','var')
    columns=1;
end;
numel=rows*columns;
indodd=1:2:2*numel;
indeven=2:2:2*numel;
nLines = 0;
while (fgets(fid) ~= -1),
  nLines = nLines+1;
end
nLines=nLines-imag(skip);
fclose(fid);
fprintf('Reading %u values from file %s.\n', nLines-1, flnm)
fid = fopen(flnm);
r=zeros(rows*columns,nLines-1);
if fsize>0
    % skip header
if skip>0
    s1 = fscanf(fid, '%s\n', real(skip));
end;
    n=0;
    while ~feof(fid)
    r0=fscanf(fid, ['(' '%f' ',' '%f' ')' '\n'], 2*numel);
    n=n+1;
    if n< nLines
    r(:,n)=r0(indodd)+1i*r0(indeven);
    end;
    end;
end;
fclose(fid);