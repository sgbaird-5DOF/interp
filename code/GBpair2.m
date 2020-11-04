function [o3_out,omega3,omega3_GBdist,varargout] = GBpair2(o1,o2,o3,varargin)
% GBPAIR2  deprecated
prec = 6;
tol = 1e-6;

method = 2;
switch method
	case 1
		%both with respect to o3
		
		%calculate distances
		[omega12,oct_sym12,zeta12,wveclist12,octonion_pair_sym_list12] = GBdist2([o1 o3],32,false);
		[omega13,oct_sym13,zeta13,wveclist13,octonion_pair_sym_list13] = GBdist2([o2 o3],32,false);
		
		wveclist3 = zeros(size(wveclist12,1),size(wveclist13,1));
		for i = 1:size(wveclist12)
			for j = 1:size(wveclist13)
				wveclist3(i,j) = wveclist12(i)+wveclist13(2);
			end
		end
		
		[omega3,minID3] = min(wveclist3);
		%find all symmetrized octonions with same omega
		minIDs12 = find(ismembertol(wveclist3,omega3,1e-6,'DataScale',1));
		
	case 2
		%both with respect to o1
		
		%calculate distances
		[omega12,oct_sym12,zeta12,wveclist12,octonion_pair_sym_list12] = GBdist2([o1 o2],32,false);
		[omega13,oct_sym13,zeta13,wveclist13,octonion_pair_sym_list13] = GBdist2([o1 o3],32,false);
end

%find all symmetrized octonions with same omega
minIDs12 = find(ismembertol(wveclist12,omega12,1e-6,'DataScale',1));
octpairsymlist12 = octonion_pair_sym_list12(minIDs12,:);

minIDs13 = find(ismembertol(wveclist13,omega13,1e-6,'DataScale',1));
octpairsymlist13 = octonion_pair_sym_list13(minIDs13,:);

%remove duplicate rows (low tol OK b.c. matching against 16 numbers)
[~,minIA] = uniquetol(round(octpairsymlist12,12),1e-3,'ByRows',true,'DataScale',1);
[~,minIA2] = uniquetol(round(octpairsymlist13,12),1e-3,'ByRows',true,'DataScale',1);

min12 = octpairsymlist12(minIA,:);
min13 = octpairsymlist13(minIA2,:);

%get omega values of combinations of second octonions from sets
%initialize
wveclist23 = zeros(1,size(min12,1)*size(min13,1));
k = 0;
for i = 1:size(min12,1)
	for j = 1:size(min13,1)
		k = k +1;
		wveclist23(k) = get_omega(min12(i,9:16),min13(j,9:16));
		min12list(k,:) = min12(i,9:16);
		min13list(k,:) = min13(j,9:16);
	end
end

%get minimum omega value within precision
[mymin,~] = min(round(wveclist23,prec));

%get corresponding octonions
myminIDs = find(abs(round(wveclist23 - mymin,prec)) < tol);
o12 = min12list(myminIDs,:); %output
o13 = min13list(myminIDs,:); %output

o12 = uniquetol(round(o12,prec),tol,'ByRows',true);
o13 = uniquetol(round(o13,prec),tol,'ByRows',true);

%arbitrarily take first octonion
if size(o12,1) > 2
	disp('')
end

% qA_0 > qB_0 convention added based on discussion with Toby Francis
if (o12(1,1) > o12(1,5)) || (size(o12,1) == 1)
	o12_out = o12(1,:);
else
	o12_out = o12(2,:);
end

if (o13(1,1) > o13(1,5)) || (size(o13,1) == 1)
	o13_out = o13(1,:);
else
	o13_out = o13(2,:);
end

%calculate distance again using GBdist (for comparison)
[omega23,oct_sym23,zeta23] = GBdist([o12_out o13_out],32,false);

%%display results
mat = [omega12;omega13;mymin;omega23];
T = array2table(mat,'VariableNames',{'values'},...
	'RowName',{'Omega12','Omega13','Omega23_pair','Omega23_GBdist'});
% 	disp([name1 '-->' name2 ', ' name1 '-->' name3])
disp(T);
end