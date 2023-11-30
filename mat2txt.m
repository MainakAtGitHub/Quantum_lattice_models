function m=mat2txt(matfile_mat2txt)
try
    load(matfile_mat2txt);
    variable_mat2txt=who;
    ismat_mat2txt=true;
catch
    load(matfile_mat2txt,'-mat');
    variable_mat2txt=who;
    ismat_mat2txt=false;
end;

%open file with write permission
if ismat_mat2txt
    fid_mat2txt = fopen([matfile_mat2txt(1:length(matfile_mat2txt)-4),'.txt'], 'w');
else
    fid_mat2txt = fopen([matfile_mat2txt,'.txt'], 'w');
end;
for i_fid_mat2txt=1:length(variable_mat2txt)
    %write a line of text
    if ~strcmp(variable_mat2txt{i_fid_mat2txt},'matfile_mat2txt')
    if isa(eval(variable_mat2txt{i_fid_mat2txt}),'char')
        fprintf(fid_mat2txt, '%s\n', [variable_mat2txt{i_fid_mat2txt},'=',eval(variable_mat2txt{i_fid_mat2txt})]);
    end;
    if isa(eval(variable_mat2txt{i_fid_mat2txt}),'float')
        fprintf(fid_mat2txt, '%s\n', [variable_mat2txt{i_fid_mat2txt},'=',num2str(eval(variable_mat2txt{i_fid_mat2txt}))]);
    end;
    end;
    %close file
end;
fclose(fid_mat2txt);
m=1;