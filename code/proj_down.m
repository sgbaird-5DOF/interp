function [projpts,usv,zeropt] = proj_down(pts,tol,usv,NV)
arguments
	pts double {mustBeFinite,mustBeReal}
	tol(1,1) double {mustBeFinite,mustBeReal} = 1e-5
	usv struct = struct.empty
	NV.nforce double = 1
	NV.nforceQ(1,1) logical = false
	NV.zeroQ(1,1) logical = true
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
% References:
%  https://www.mathworks.com/matlabcentral/answers/352830
%
%--------------------------------------------------------------------------
% unpackage name value pairs
nforce = NV.nforce;
nforceQ = NV.nforceQ;

% --if zeropt is requested as output, then zeroQ == true
if nargout == 3
	zeroQ = true;
else
	zeroQ = NV.zeroQ;
end

%dimensionality
d = size(pts,2);

if nforce >= d
	error(['nforce should be less than d == ' int2str(size(pts,2))])
end

if ~isempty(usv)
	%unpackage usv
	V = usv.V;
	S = usv.S;
	avg = usv.avg;
	
	%projection
	projpts = (pts-avg)/V';
	
	%number of degenerate dimensions
	if nforceQ
		ndegdim = nforce;

	elseif size(S,1) == size(S,2)
		ndegdim = sum(abs(diag(S)) < tol);
		
	else
		%check for columns of all zeros within tolerance
		ndegdim = sum(all(S < tol));
	end
	
	if all(abs(projpts(:,end-ndegdim+1:end)) < tol,'all')
		%remove last column(s)
		projpts = projpts(:,1:end-ndegdim);
		
	elseif nforceQ
		projpts = projpts(:,1:end-ndegdim);
		disp(['Nonzero last column. E.g. ' num2str(pts([1 2],end)) ...
			'. Forcing projection ' int2str(ndegdim) ' dimensions.'])
		
	elseif ~nforceQ
		projpts = pts;
		usv = struct.empty;
		%not sure if I should have a constant, non-zero last column be OK
		if size(pts,1) > 3
			n = 3;
		else
			n = 1;
		end
		warning(['Nonzero last column. E.g. ' num2str(pts(1:n,end).') '. Setting projpts == pts'])
	end
	
elseif isempty(usv)
	% make a non-empty struct with no fields
	usv(1) = struct();
	
	%take average of points
	avg = mean(pts);
	
	%project to d-1 dimensional space
	[U,S,V] = svd(pts-avg,0);
	
	usv.U = U;
	usv.S = S;
	usv.V = V;
	usv.avg = avg;
	
	%number of degenerate dimensions
	if nforceQ
		ndegdim = nforce;
		
	elseif size(S,1) == size(S,2)
		ndegdim = sum(abs(diag(S)) < tol);
		
	else
		%check for columns of all zeros within tolerance
		ndegdim = sum(all(S < tol));
	end
	
	if (ndegdim > 0) || nforceQ
		%project to lower dimension (i.e. rotation)
		projpts = U*S(:,1:d-ndegdim);
		if (ndegdim == 0) && nforceQ
			warning(['ndegdim == 0, tol == ' num2str(tol) ...
				', min(diag(S)) == ' num2str(min(diag(S))) ...
				', max(diag(S)) == ' num2str(max(diag(S))) ...
				'. Forcing projection ' int2str(ndegdim) ' dimensions'])
		end
	else
		warning(['ndegdim == 0, tol == ' num2str(tol) ...
			', min(diag(S)) == ' num2str(min(diag(S))) ...
			', max(diag(S)) == ' num2str(max(diag(S))) ...
			'. Setting projpts == pts'])
		projpts = pts;
		usv = struct.empty;
	end
end

if zeroQ
	zeropt = (zeros(1,d)-avg)/V';
	zeropt = zeropt(1:d-ndegdim);
	projpts = projpts-zeropt;
	
	usv.zeropt = zeropt;
end

end %proj_down

%-----------------------------CODE GRAVEYARD-------------------------------
%{
	% ndegdim = sum(abs(diag(S)) < tol);

% if isempty(nforce)
% 	nforce = 1;
% 	nforceQ = false;
% else
% 	nforceQ = true;
% end

% 	[U,S,V]=svd(bsxfun(@minus,pts,avg),0);


%}