function [projpts,usv] = proj_down(pts,tol,usv,nforce)
arguments
	pts double {mustBeFinite,mustBeReal}
	tol(1,1) double {mustBeFinite,mustBeReal} = 1e-6
	usv struct = struct.empty
	nforce(1,1) double {mustBeInteger} = 1
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Usage:
%		[projpts,usv] = proj_down(pts);
%
%		[projpts,usv] = proj_down(pts,usv);
%
% Date:
%
% Description: project down by removing null dimensions (i.e. via rotation)
%
% Inputs:
%
% Outputs:
%
% Dependencies:
%
% Notes:
%	Small numerical errors can accumulate by calling proj_down & proj_up
%	repeatedly with different usv matrices. Shouldn't be an issue if you
%	keep use the same usv matrices.
%
%--------------------------------------------------------------------------
%dimensionality
d = size(pts,2);



if nforce >= d
	error('nforce should be less than d == size(pts,2)')
elseif nforce > 1
	nforceQ = true;
end

if ~isempty(usv)
	%unpackage
	V = usv.V;
	avg = usv.avg;
	
	%projection
	projpts = (pts-avg)/V';
	
	if all(abs(projpts(:,end-nforce+1:end)) < tol,'all') || nforceQ
		%remove last column
		projpts = projpts(:,1:end-nforce);
		
	elseif ~nforceQ
		projpts = pts;
		%not sure if I should have a constant, non-zero last column be OK
		warning(['Nonzero last column. E.g. ' num2str(pts([1 2],end)) '. Setting projpts == pts'])
	end
	
else
	%take average of points
	avg = mean(pts);
	
	%project to d-1 dimensional space
	[U,S,V]=svd(bsxfun(@minus,pts,avg),0);
	usv.U = U;
	usv.S = S;
	usv.V = V;
	usv.avg = avg;
	
	%number of degenerate dimensions
	ndegdim = sum(abs(diag(S)) < tol);
	
	if (ndegdim > 0) || nforceQ
	%project to lower dimension (i.e. rotation)
		projpts = U*S(:,1:d-nforce);
		if (ndegdim == 0) && nforceQ
			warning(['ndegdim == 0, tol == ' num2str(tol) ...
				', min(diag(S)) == ' num2str(min(diag(S))) ...
				', max(diag(S)) == ' num2str(max(diag(S))) ...
				'. Forcing projection ' int2str(nforce) ' dimensions'])
		end
	else
		warning(['ndegdim == 0, tol == ' num2str(tol) ...
			', min(diag(S)) == ' num2str(min(diag(S))) ...
			', max(diag(S)) == ' num2str(max(diag(S))) ...
			'. Setting projpts == pts'])
		projpts = pts;
	end
end

end %proj_down

%-----------------------------CODE GRAVEYARD-------------------------------
%{
	% ndegdim = sum(abs(diag(S)) < tol);
%}