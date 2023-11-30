function [beginstring, number,endstring]=getnumber(filename,initialstring)
% split filename after initial string in numerical value and the remaining
% string
k=strfind(filename,initialstring);
endstring=filename(k+length(initialstring):length(filename));
numericcharacters='1234567890-+.eE';
numberstring_new=[];
while sum(numericcharacters==endstring(1))
    numberstring_new=[numberstring_new,endstring(1)];
    endstring=endstring(2:length(endstring));
    if isempty(endstring)
        break
    end;
end;
number=str2double(num2str(numberstring_new));
if isempty(number)
    zpos=input(['could not find correct numeric value after ', initialstring,' in string ',filename,' , please enter: '],'s');    
end;
beginstring=filename(1:k-1);