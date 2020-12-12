%import data
addpath(genpath('.'))
fname = 'Supplementary data - GB-Energy for TIGB&TWGB in Part II.txt';
T = fileread(fname);
% T = {horzcat(T{:,:})}
% T = fillmissing(T,'constant',[])
%split into header/data
nheaders = 16;
headerinfo = T(1:nheaders);
data = T(nheaders+1:end);
%unpack data
misaxtxt = data(1:5:end);
misangletxt = data(2:5:end);
csltxt = data(3:5:end);
tiltGBEtxt = data(4:5:end);
twistGBEtxt = data(5:5:end);



%% CODE GRAVEYARD
%{
% fid = fopen(fname)
% i = 0;
% while ~feof(fid)
%     i = i+1;
%     a = fgetl(fid);
% end
% T = {};
% fclose(fid);
% fid = fopen(fname);
% while ~feof(fid)
%     T = [T,fgetl(fid)];
% end
%}