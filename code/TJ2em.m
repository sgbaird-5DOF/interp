function [e1,e2,e3,m1,m2,m3,intIDs] = TJ2em(EAs,norms)
arguments
    EAs(:,3,3) double {mustBeReal,mustBeFinite}
    norms(:,3,3) double {mustBeReal,mustBeFinite}
end
% TJ2EM  triple junctions to sorted euler (e1,e2,e3) and BP normals (m1,m2,m3) and interleaf IDs (intIDs)

%number of
%--triple junctions
nTJ = size(EAs,1);
%--grain boundaries
nGB = nTJ*3;

%IDs for interleaving vectors
%e.g. X = [xi_1^x xi_1^y xi_1^z xi_2^x xi_2^y xi_2^z ... ], see [1]
id1 = 1:3:nGB;
id2 = 2:3:nGB;
id3 = 3:3:nGB;
intIDs = [id1 id2 id3]; %interleaf IDs, i.e. [1 4 7 ... 2 5 8 .. 3 6 9 ... ]

% catenate and squeeze helper function (apply to each grain and normal)
catsqz = @(mat,i) squeeze(vertcat(mat(:,i,:)));

%extract euler angles
e1 = catsqz(EAs,1); % i.e. e1 = squeeze(vertcat(EAs(:,1,:));
e2 = catsqz(EAs,2);
e3 = catsqz(EAs,3);

%extract BP normals
m1 = catsqz(norms,1); % i.e. m1 = squeeze(vertcat(norms(:,1,:));
m2 = catsqz(norms,2);
m3 = catsqz(norms,3);
end