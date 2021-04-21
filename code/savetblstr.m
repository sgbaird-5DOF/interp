function tblstr = savetblstr(tbl,savename,savefolder,caption,nv)
arguments
    tbl
    savename = 'tblstr-test'
    savefolder = '.'
    caption = '...'
    nv.printQ(1,1) logical = true
end
%SAVETBLSTR  compile a \begin{table}...\end{table} string and save to .tex file
printQ = nv.printQ;
% tblstr = latexTable(struct('data',tbl,'tableLabel',savename,'tableCaption',caption,'booktabs',true));
tblstr = [...
    '\\begin{table}[%s]\n'...
    '\t\\centering\n'...
    '\t\\label{%s}\n'...
    '\t\\caption{%s}\n'...
    
    '\\end{table}'];
tblstr = strjoin(tblstr,'\n');
savepath = fullfile(savefolder,[savename,'.tex']);
f = fopen(savepath,'w');
fprintf(f,'%s',tblstr);
tblstr = fileread(savepath);
if printQ
    disp(tblstr);
end
fclose(f);
end