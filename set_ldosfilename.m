if ~(exist('sqstring','var'))
    sqstring='';
end;
LDOSfileName0 = [casestring,'_Vimp_', num2str(Vimp),  '_N_', num2str(N)];
LDOSfileName = [LDOSfileName0 , '_M_', num2str(M),'_ita_', num2str(ita),'_e_',num2str(E),sqstring];