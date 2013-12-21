% script to read the file read_input_file and load all variables that are
% given there (with support of double variables and strings)
fid = fopen(read_input_file);
tline = fgetl(fid);
while ischar(tline)
    % display the text (for debugging)
    % disp(tline)
    % manipulate the string
    % 1) remove part after a % sign (to write comments)
    k_position = strfind(tline, '%');
    if ~(isempty(k_position))
        tline=tline(1:k_position-1);
    end;
    % 2) remove all spaces
     tline(isspace(tline)) = [];
    % possibly also remove other characters?
    %  stringname(ismember(stringname,' ,.:;!')) = [];
    % 3) search for pattern <variablename>=<value>
    k_position= strfind(tline, '=');
    if ~(isempty(k_position))
        variablename=tline(1:k_position-1);
        if numel(k_position)==1
            value=tline(k_position+1:length(tline));
        else
            value=tline(k_position(1)+1,k_position(2));
        end
    number=str2num(value);
    if ~isempty(number)
        % 4) a) <value> represents a number
        disp([variablename,'=',num2str(number)])
        eval([variablename,'=number;']);
    else
        %    b) <value> represents a string
        disp([variablename,'=',value])
        eval([variablename,'=value;']);
    end
    end;
    tline = fgetl(fid);
end
fclose(fid);
% clean up
clear value variablename tline k_position number fid ans