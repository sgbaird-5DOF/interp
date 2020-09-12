%parseMyFiles_test
[fpathshort,nameExt,plist] = parseReqFiles('interp5DOF.m');
fpathshort{1}
nameExt{1}
plist

T = table(['---'; nameExt],['---'; fpathshort],'VariableNames',{'name','folder'});
disp(T)
writetable(T,'interp5DOF_requiredFiles.txt','Delimiter','|')

T2 = table(strcat('1. [',nameExt,'](',fpathshort,nameExt,')'),'VariableNames',{'### names'});
disp(T2)
writetable(T2,'interp5DOF_requiredFilesGitHub.txt','Delimiter','|','WriteVariableNames',false)