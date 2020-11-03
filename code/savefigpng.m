function savefigpng(folder,fname)

fpath = fullfile(folder,fname);
savefig(fpath)
print(fpath,'-dpng','-r300')

end