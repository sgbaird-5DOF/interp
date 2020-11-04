function write_video(images,movname)
% WRITE_VIDEO  write a set of images to a video named movname with some defaults
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