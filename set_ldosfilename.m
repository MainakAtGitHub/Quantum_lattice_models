if ~(exist('sqstring','var'))
    sqstring='';
end;
if ischar(xrange)
        sqstring=[sqstring,xrange];
end;
if ~(exist('en','var'))
    en=1;
end;
LDOSfileName0 = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N)];
LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_ita_', num2str(ita),'_e_',num2str(E(en)),sqstring];