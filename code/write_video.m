function write_video(images,movname)

%close all %not sure if you need this

mymovie = VideoWriter(movname); %,'MPEG-4');
mymovie.FrameRate = 1;
mymovie.Quality = 95;
open(mymovie)
writeVideo(mymovie,[images{1}]);
for i = 2:length(images)
	writeVideo(mymovie,[images{i}]);
end
close(mymovie)

end