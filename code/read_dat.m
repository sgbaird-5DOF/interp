function [TJs,EAs,norms] = read_dat(fpath,nheaderlines)
arguments
    fpath char = '../../TJ2GBE/TJdata/triples_30000.dat'
    nheaderlines(1,1) double = 0
end
% READ_DAT  extract TJs, EAs, and norms from triples.dat file. (in-progress, 2020-12-05)
% Mirrors https://github.com/Yufeng-shen/TJ2GBE/blob/master/Src/Python/Reconstruction.py
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-12-05
%
%     Input: triples.dat, wrote from the fortran program Torq_gen
%                 size=[numTJ*8,]
%                 In each group, the data is [TJ directon, EA1, GB1, EA2, GB2, EA3, GB3]
%     Output: TJs, direction of the triple junctions
%                 size = [numTJ, 3]
%             EAs, the EA angles of the 3 grains at a TJ
%                 size = [numTJ, 3, 3]
%             norms, normal direction of the 3 GB at a TJ
%                 size = [numTJ, 3, 3]
%
% 1 - TJnum             (lines 1, 8, etc.)
%   2 3 4 - TJs         (lines 2, 9, etc.)
%   5 6 7 - EAs         (lines 3, 10, etc.)
%   8 9 10 - EAs        (lines 4, 11, etc.)
%   11 12 13 - EAs      (lines 5, 12, etc.)
%   14 15 16 - norms    (lines 6, 13, etc.)
%   17 18 19 - norms    (lines 7, 14, etc.)
%   20 21 22 - norms    (lines 8, 15, etc.)
%--------------------------------------------------------------------------
data = importdata(fpath,' ',nheaderlines);
TJs = [data(2:22:end),data(3:22:end),data(4:22:end)]; %lines 2, 9, etc.
EAs = [...
    data(5:22:end),data(6:22:end),data(7:22:end); %lines 3, 10, etc.
    data(8:22:end),data(9:22:end),data(10:22:end); %lines 4, 11, etc.
    data(11:22:end),data(12:22:end),data(13:22:end)]; %lines 5, 12, etc.

norms = [...
    data(14:22:end),data(15:22:end),data(16:22:end); %lines 6, 13, etc.
    data(17:22:end),data(18:22:end),data(19:22:end); %lines 7, 14, etc.
    data(20:22:end),data(21:22:end),data(22:22:end)]; %lines 8, 15, etc.
% unecessary, but useful as a double check
TJnum = data(1:22:end); %lines 1, 8, etc.
nTJs = length(TJnum); %#ok<NASGU>

end
%% CODE GRAVEYARD
%{
% nTJs = length(TJnum);
% TJ = zeros(nTJs,3);
%}