%parseMyFiles_test
fname = 'interp5DOF.m';
f = dir(fname);
folder = f(1).folder;
addpath(genpath(folder));
sep = '/';
[fpathshort,nameExt,plist] = parseReqFiles(fname,'interp*','sep',sep);
fpathshort{1}
nameExt{1}
plist

T = table(['---'; nameExt],['---'; fpathshort],'VariableNames',{'name','folder'});
disp(T)
writetable(T,'interp5DOF_requiredFiles.txt','Delimiter','|')

T2 = table(strcat('1. [',nameExt,'](',fpathshort,nameExt,')'),'VariableNames',{'### names'});
disp(T2)
writetable(T2,'interp5DOF_requiredFilesGitHub.txt','Delimiter','|','WriteVariableNames',false)
