function [ax,D,nnpts,idx] = nnhist(pts,dtype)
arguments
	pts
	dtype char {mustBeMember(dtype,{'omega','norm','alen'})} = 'omega'
end
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-07-27
%
% Description: nearest neighbor distance histogram
% 
% Inputs:
%  pts - rows of points (euclidean)
%
%  dtype - distance type
%
% Outputs:
%  creates a histogram figure
%
% Usage:
%	pts = sqrt2norm(normr(rand(1000,8))); %octonion case (omega), input must have norm == sqrt(2)
%  nnhist(pts)
%
%  pts = nnhist(rand(1000,5)); % euclidean distance
%  nnhist(pts)
%
% Dependencies:
%  get_omega.m
%
% Notes:
%  *
%--------------------------------------------------------------------------
%get nearest neighbor IDs and euclidean distances
[idxtmp,Dtmp] = knnsearch(pts,pts,'K',2);
idx = idxtmp(:,2);

%nearest neighbor pts
nnpts = pts(idx,:);

figure

%get distances for plotting
switch dtype
	case 'omega'
		Drad = get_omega(pts,nnpts);
		D = rad2deg(Drad);
		xlbl = 'NN \omega (deg)';
	case 'norm'
		D = Dtmp(:,2);
		xlbl = 'NN euclidean octonion distance';
	case 'alen' %general arc length formula
		Drad = real(acos(dot(pts,nnpts,2)));
		D = rad2deg(Drad);
		xlbl = 'NN \omega = cos^{-1}(p_1 \cdot p_2) (deg)';
end

%plotting
ax = histogram(D);
xlabel(xlbl)
ylabel('counts')

%---------------------------------CODE GRAVEYARD---------------------------
%{
% pts = uniquetol(pts,'ByRows',true);

%}