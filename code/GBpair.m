function [o3_out,omega3,omega3_GBdist,varargout] = GBpair(o1,o2,o3,method)
arguments
	o1(:,8) double {mustBeFinite,mustBeReal}
	o2(:,8) double {mustBeFinite,mustBeReal}
	o3 %can be empty if method == standard
	method string {mustBeMember(method,'standard','pairwise')} = 'standard'
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-17
%
% Description: Method 1: 
%						Find o3 that has the minimum summed distances, o1-->o3,
%					and o2-->o3.
%					
%					Method 2:
%						Find the symmetrically equivalent minimum o2's and
%					o3's each with respect to o1 and take the (o2,o3) pair from
%					these that have a minimum with respect to each other (i.e.
%					o2-->o3)
%
% Inputs:
%		o1			===	first octonion, stays constant no matter which method
%							is used
%
%		o2			===	second octonion, stays constant in method 1, allowed
%							to change in method 2
%
%		o3			===	third octonion, always gets symmetrized
%
%		method	===	which method to use (see description)
%
% Outputs:
%
%		o3_out	===	symmetrized third octonion
%
%		omega3	===	summed omega distances of o1-->o3_out and o2-->o3_out
%							(method 1). Alternatively, distance from
%							o2_out-->o3_out
%
%		omega3_GBdist	===	Summed distances of two standard GBdist()
%									function calls, o1-->o3_out and o2-->o3_out
%									(method 1). Alternatively, standard GBdist()
%									function call using o2_out and o3_out (method
%									2). Can be considered as a "ground truth" to
%									compare against to omega3.
%
%		o2_out	===	symmetrized second (optional output, method 2 only)
%
% Dependencies:
%
%--------------------------------------------------------------------------
if strcmp(method,'pairwise') && isempty(o3)
	error("o3 must be specified if pairwise method is used. Did you mean GBpair(o1,o2)?")
end

prec = 6;
tol = 1e-6;

switch method
	case 'standard'
		% o3 with respect to o1 and o2, output o3
		
		%calculate distances
		[omega1,oct_sym1,zeta1,wveclist1,octonion_pair_sym_list1] = GBdistEucl([o1 o3],32,false);
		[omega2,oct_sym2,zeta2,wveclist2,octonion_pair_sym_list2] = GBdistEucl([o2 o3],32,false);
		
	case 'pairwise'
		%both with respect to o1, output o2 and o3
		
		%calculate distances
		[omega1,oct_sym1,zeta1,wveclist1,octonion_pair_sym_list1] = GBdistEucl([o1 o2],32,false);
		[omega2,oct_sym2,zeta2,wveclist2,octonion_pair_sym_list2] = GBdistEucl([o1 o3],32,false);
end

%find all symmetrized octonions with same minimum omega
minIDs1 = find(ismembertol(wveclist1,omega1,tol,'DataScale',1));
octpairsymlist1 = octonion_pair_sym_list1(minIDs1,:);
% wveclist1 = wveclist1(minIDs);

minIDs2 = find(ismembertol(wveclist2,omega2,tol,'DataScale',1));
octpairsymlist2 = octonion_pair_sym_list2(minIDs2,:);
% wveclist2 = wveclist2(minIDs);

%remove duplicate rows (low tol OK b.c. matching against 16 numbers)
[~,minIA1] = uniquetol(round(octpairsymlist1(:,9:16),prec),tol,'ByRows',true,'DataScale',1);
[~,minIA2] = uniquetol(round(octpairsymlist2(:,9:16),prec),tol,'ByRows',true,'DataScale',1);

%extract unique rows
min1 = octpairsymlist1(minIA1,:);
min2 = octpairsymlist2(minIA2,:);

switch method
	case 'standard'
		skipQ = true;
		if skipQ
			omega3 = omega2;
			o3_out = octpairsymlist2(1,9:16);
			
			%calculate distance again using GBdist (for comparison)
			omega3_GBdist = omega2;
			return
		end
		%output o3. o1 and o2 kept constant
		
		%get omega values corresponding to unique rows
% 		w1 = wveclist1(minIA1);
% 		w2 = wveclist2(minIA2);
		
		%find o3 values that are the same
		[min3,ia,ib] = intersect(min1(:,9:16),min2(:,9:16),'rows');
		
		%take o3's with minimum summed distance
% 		wveclist3 = w1(ia)+w2(ib);
		
% 		[omega3,ic] = min(wveclist3);
% 		minIDs3 = find(ismembertol(wveclist3,omega3,tol,'DataScale',1));
% 		symlist3 = min3(minIDs3,:);
		
		symlist3 = min3;
		
		if size(symlist3,1) > 2
			warning(['size(symlist3,1) == ' int2str(size(symlist3,1))])
		end
		
		% parse output
		% qA_0 > qB_0 convention added based on discussion with Toby Francis
		if (symlist3(1,1) > symlist3(1,5)) || (size(symlist3,1) == 1)
			o3_out = symlist3(1,:);
		else
			o3_out = symlist3(2,:);
		end
		
		%calculate distance again using GBdist (for comparison)
		wTemp1 = GBdist([o1 o3_out],32,false,false);
		wTemp2 = GBdist([o2 o3_out],32,false,false);
		omega3_GBdist = wTemp1+wTemp2;
		
		%display results
% 		mat = [omega1;omega2;omega3];
% 		T = array2table(mat,'VariableNames',{'values'},...
% 			'RowName',{'Omega1','Omega2','Omega3'});
% 		disp(T);
		
	case 'pairwise'
		%output o2 and o3 with respect to o1
		
		%initialize
		npts3 = size(min1,1)*size(min2,1);
		wveclist3 = zeros(1,npts3);
		minlist1 = zeros(npts3,8);
		minlist2 = minlist1;
		k = 0;
		for i = 1:size(min1,1)
			for j = 1:size(min2,1)
				%get omega values of combinations of second octonions
				k = k+1;
				wveclist3(k) = get_omega(min1(i,9:16),min2(j,9:16));
				minlist1(k,:) = min1(i,9:16);
				minlist2(k,:) = min2(j,9:16);
			end
		end
		
		%get minimum omega value within precision
		omega3 = min(round(wveclist3,prec));
		
		%get corresponding octonions
		ids = find(abs(round(wveclist3-omega3,prec)) < tol);
		
		o12 = minlist1(ids,:);
		o13 = minlist2(ids,:);
		
		o12 = uniquetol(round(o12,prec),tol,'ByRows',true,'DataScale',1);
		o13 = uniquetol(round(o13,prec),tol,'ByRows',true,'DataScale',1);
		
		%arbitrarily take first octonion
		if size(o12,1) > 2
			disp('')
		end
		
% 		% qA_0 > qB_0 convention added based on discussion with Toby Francis
% 		if (o12(1,1) > o12(1,5)) || (size(o12,1) == 1)
% 			o2_out = o12(1,:);
% 		else
% 			o2_out = o12(2,:);
% 		end
% 		
% 		if (o13(1,1) > o13(1,5)) || (size(o13,1) == 1)
% 			o3_out = o13(1,:);
% 		else
% 			o3_out = o13(2,:);
% 		end
		
		o2_out = o12(1,:);
		o3_out = o13(1,:);
		
		%calculate distance again using GBdist (for comparison)
		[omega3_GBdist,oct_sym3,zeta3] = GBdist([o2_out o3_out],32);
		
		%package output
		varargout{1} = o2_out;
		varargout{2} = omega3_GBdist;
		
		%display results
% 		mat = [omega1;omega2;omega3;omega3_GBdist];
% 		T = array2table(mat,'VariableNames',{'values'},...
% 			'RowName',{'Omega1','Omega2','Omega3','Omega3_GBdist'});
% 		disp(T);
end

if (method == 2) && exist('o2_out','var') == 0
	warning('o2 and/or omega3_GBdist output specified, but method 1 selected. Remove output or change to method 2.')
end


end %GBpair
