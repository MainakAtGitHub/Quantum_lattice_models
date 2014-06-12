function str=my_int2str(in)
% code usable in conjunction with code generator
str='';
minus='';
if in <0
    minus='-';
    in=-in;
end;
while in >0
   str=[int2str_single(mod(in,10)),str];
   % note floor and fix don't work with the code generator!
   in=floor(in/10);
end;
if isempty(str)
    str='0';
end;
str=[minus,str];
function c=int2str_single(in)
switch in
    case 0
        c='0';
    case 1
        c='1';
    case 2
        c='2';
    case 3
        c='3';
    case 4
        c='4';
    case 5
        c='5';
    case 6
        c='6';
    case 7
        c='7';
    case 8
        c='8';
    otherwise
        c='9';
end
