function figstr = savefigstr(caption,savename,savefolder,figopt,nv)
arguments
    caption = get_captionlist()
    savename = 'get_figstr-test'
    savefolder = '.'
    figopt = 'scale=1'
    nv.printQ(1,1) logical = true
end
%SAVEFIGSTR  compile a \begin{figure}...\end{figure} string and save to .tex file
printQ = nv.printQ;
figstr = [...
    '\\begin{figure}[h]\n'...
    '\t\\centering\n'...
    '\t\\includegraphics{%s}[%s]\n'...
    '\t\\label{%s}\n'...
    '\t\\caption{%s}\n'...
    '\\end{figure}'];
savepath = fullfile(savefolder,[savename,'.tex']);
f = fopen(savepath,'w');
fprintf(f,figstr,[savename '.png'],figopt,['fig:' savename],caption);
figstr = fileread(savepath);
if printQ
    disp(figstr);
end
fclose(f);
end