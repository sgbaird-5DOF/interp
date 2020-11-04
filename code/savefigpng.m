function savefigpng(folder,fname)
% SAVEFIGPNG  save a figure and print the figure as 300 DPI .png
fpath = fullfile(folder,fname);
savefig(fpath)
print(fpath,'-dpng','-r300')

end