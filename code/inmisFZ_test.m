%inmisFZ test
clear; close all

S = load('misFZfeatures.mat');

%catenate fields from dlist
vnames = fields(S.dlist);
nvars = length(vnames);
dlist = zeros(nvars,3);
for i = 1:length(vnames)
	vname = vnames{i};
	dlist(i,:) = S.dlist.(vname);
end

%get constraints & call inmisFZ
[A,b] = misFZcon();
tol = 1e-12;
TF = inmisFZ(dlist,A,b,tol);

if all(TF)
	disp('all misFZfeatures were considered inside of misFZ')
else
	disp(['# misFZfeatures outside of misFZ: ' int2str(sum(~TF)) ', within tol: ' num2str(tol)])
end



%---------------------------------CODE GRAVEYARD---------------------------
%{
geomchecktypeList = {...
	'OAB','OBCE','OADE','CDE',... surfaces
	'OB','CE','ED','OE','OA','AC',... lines
	'B','E','A','C','O'}; %points
%}