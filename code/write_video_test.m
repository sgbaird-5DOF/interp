clear; close all

vidQ = true;

savefolder = '.';

%generate plots and save pictures
for i = 1:3
	plot(i:i^2+1);
	fnames{i} = fullfile(savefolder,['plot',int2str(i) '.png']);
	print('-r150',fnames{i},'-dpng');
	images{i} = imread(fnames{i});
end

if vidQ
	movname = fullfile(savefolder,'movie_name.avi');
	write_video(images,movname)
end

%Mmovie = immovie([final_images{:}]);
%implay(mymovie)