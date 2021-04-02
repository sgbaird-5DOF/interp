function savefigpng(folder,fname,cropPos,manualcrop)
arguments
    folder char
    fname char
    cropPos double = []
    manualcrop(1,1) logical = false
end
% SAVEFIGPNG  save a figure and print the figure as 300 DPI .png
fpath = fullfile(folder,fname);
savefig(fpath)
% print(fpath,'-dpng','-r300')
exportgraphics(gcf,[fpath '.png'],'Resolution',300)
if ~isempty(cropPos) || manualcrop
    pngpath = [fpath,'.png'];
    I = imread(pngpath);
    if manualcrop
        paperfigure()
        imshow(I)
        Icropped = imcrop;
    else
        Icropped = imcrop(I,cropPos);
    end
    imwrite(Icropped,pngpath,'XResolution',11811,'YResolution',11811,'ResolutionUnit','meter')
end

%% CODE GRAVEYARD
%{
%     paperfigure(NV.nrows,NV.ncols);

%     imshow(Icropped,'InitialMagnification',1)
%     fig2 = gcf;
%     fig2.Position = pos;
%     print(fpath,'-dpng','-r300')

    NV.nrows(1,1) double = 1
    NV.ncols(1,1) double = 1

%}