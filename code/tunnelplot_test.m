%%
mdlnum = 1;
mdl = mdlcat(mdlnum);
npts = mdl.mesh.npts;
ids = randi(npts,18,1);
paperfigure(3,3);
cstAQ = true;
for i = 1:2:length(ids)
    nexttile
    if cstAQ
%         A = mdl.mesh.pts(ids(1),:);
        A = mdl.oref;
    else
        A = mdl.mesh.pts(ids(i),:);
    end
    B = mdl.mesh.pts(ids(i+1),:);
    tunnelplot(mdl,A,B);
end