function [five_out,olist_out] = tofiveFZ(five,varargin)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-21
%
% Description: Take 5DOF data and rotate it into misorientation and
% boundary plane fundamental zones.
%
% Inputs:
%		five - struct array of 5DOF parameters
%
% Outputs:
%		five_out - struct array of 5DOF parameters that have been converted
%		to FZ representation
%
% Dependencies:
%
% Notes:
%		Assumes that eadch spherical triangle arc falls only along longitudes
%		or latitudes (i.e. not a mix of both), which is why there is a rotation
%--------------------------------------------------------------------------
if nargin == 2
	olist = varargin{1};
else
	olist = GBfive2oct(vertcat(five.q),vertcat(five.nA));
end

npts = length(five);
infiveFZ = zeros(1,npts);
olist_out = zeros(npts,8);

%initialize
five_out(npts) = struct;
five_out(1).q = [];
five_out(1).d = [];
five_out(1).nA = [];
five_out(1).geometry = '';

%textwaitbar setup
waitbarQ = false;
if waitbarQ
	D = parallel.pool.DataQueue;
	afterEach(D, @nUpdateProgress);
	N=npts;
	p=1;
	reverseStr = '';
	nreps = floor(N/100);
	nreps2 = nreps;
	
end

% 	function nUpdateProgress(~)
% 		percentDone = 100*p/N;
% 		msg = sprintf('Percent done: %3.1f ', percentDone); %Don't forget this semicolon
% 		fprintf([reverseStr, msg]);
% 		reverseStr = repmat(sprintf('\b'), 1, length(msg));
% 		p = p + nreps;
% 	end

for i = 1:npts %parfor compatible
% for i = 2
	%% setup
	%unpackage point
	o = olist(i,:);
	q = five(i).q;
	d = five(i).d;
	nA = five(i).nA;
	geometry = five(i).geometry;
	
	%generate symmetrically equivalent octonions
	osyms = osymsets(o,32);
	osyms = osyms{1};
	
	%convert to 5DOF
	symfive = GBoct2five(osyms,'parforQ',false,'disQ',false); %call to disorientation is expensive
	
	[~,ia] = uniquetol(round([vertcat(symfive.q) vertcat(symfive.nA)],6),1e-3,'ByRows',true,'DataScale',1);
	symfive = symfive(ia);
	
	%compute rodrigues vectors
	% 	dlist = vertcat(symfive.d);
	qlist = vertcat(symfive.q);
	dlist = q2rod(qlist);
	
	% get rotation matrix
	nint = 5;
	ctrcuspQ = false;
	[pts,A,R,~,~,sphpts] = meshBP(q,nint,ctrcuspQ,geometry);
	%pts are in the lab frame. sphpts are in the sample frame
	vertices = vertcat(sphpts{:});
	
	symfive1 = symfive;
	
	%% find misFZ IDs
	%take only those in misFZ
	% 	misFZ_ids = InCubicFZ(dlist);
	
	misFZ_ids = false(size(qlist,1),1);
	for j = 1:size(qlist,1)
		geom = findgeometry(qlist(j,:));
		if ~strcmp(geom,'exterior')
			misFZ_ids(j) = true;
		end
	end
	
	symfive = symfive(misFZ_ids);
	osyms = osyms(misFZ_ids,:);
	
	%% find BPFZ IDs
	%test for points in BPFZ
	nAlist = vertcat(symfive.nA); %lab frame (?)
	BPFZ_ids = inBPFZ(nAlist,vertices,R); %nA list gets rotated to sample frame using R
	
	%take only first one found if multiple? how to make sure at least one is found?
	symfive = symfive(BPFZ_ids);
	osyms = osyms(BPFZ_ids,:);
	
	%% plotting
	if ismember(i,[1 2])
		plotQ = true;
	else
		plotQ = false;
	end
	S = var_names(symfive1,dlist,pts,nAlist,nA,A,R,plotQ,misFZ_ids,geometry);
	if plotQ
		myplot(S);
	end
	
	%% final output
	%take only those in BPFZ
	five_out(i) = symfive;
	olist_out(i,:) = osyms;
	misFZ_ids2 = find(misFZ_ids);
	
	%find common indices (i.e. in 5DOF FZ)
	infiveFZ(i) = misFZ_ids2(BPFZ_ids);
	
	
	%% text waitbar
	if waitbarQ
		if mod(i,nreps2) == 0
			send(D,i);
		end
	end
end

end

%%
%------------------------------HELPER FUNCTIONS----------------------------
function myplot(S)

vars = fields(S);

for i = 1:length(vars)
	var = vars{i};
	temp = S.(var); %#ok<NASGU> %temporary value of vName
	evalc([var '= temp']); %assign temp value to a short name
end

if plotQ
	clf(figure(1)); fig = figure(1);
	fig.Position = [211.5000  216.0000  876.0000  537.5000];
	tiledlayout(2,3)
	t = nexttile(1);
	%misFZ
	plotFZrodriguez_vtx();
	xlims = t.XLim;
	ylims = t.YLim;
	zlims = t.ZLim;
	hold on
end

if plotQ
	%misFZ with all symmetrically equivalent points
	t=num2cell(dlist,1);
	plot3(t{:},'.')
	K = sphconvhulln(dlist);
	
	%BPFZ lab frame triangulation of all corresponding normals
	nexttile
	nAtemp = vertcat(symfive1.nA);
	nAtemp1 = (R\nAtemp.').';
	t = num2cell(nAtemp1,1);
	K = sphconvhulln(nAtemp1);
	trisurf(K,t{:})
	title('lab frame')
	axis equal
	
	%BPFZ sample frame triangulation of all corresponding normals
	nexttile
	nAtemp2 = nAtemp;
	t = num2cell(nAtemp2,1);
	K = sphconvhulln(nAtemp2);
	trisurf(K,t{:})
	title('sample frame')
	axis equal
end

if plotQ
	%misFZ of points inside misFZ
	nexttile;
	plotFZrodriguez_vtx();
	t = num2cell(dlist(misFZ_ids,:),1);
	plot3(t{:},'*')
end

if plotQ
	pts = uniquetol(vertcat(pts.sub),'ByRows',true);
	
	%BPFZ lab frame
	%normals
	nexttile
	% 		t = (R\nAlist.').';
	t = nAlist;
	t = num2cell(1.2*t,1);
	plot3(t{:},'r*')
	
	hold on
	
	%---------true point--------
	t = num2cell(1.3*nA);
	plot3(t{:},'ko')
	
	%---------true point--------
	t = num2cell(1.3*(R\nA.').');
	plot3(t{:},'bo')
	
	%---------true point--------
	t = num2cell(1.3*(R*nA.').');
	plot3(t{:},'co')
	
	sphere(40);
	
	scl = 1.2;
	x = A(:,1)*scl;
	y = A(:,2)*scl;
	z = A(:,3)*scl;
	
	t = num2cell(A,2);
	
	w = [0,0,0];
	quiver3(w,w,w,t{:},0,'linewidth',1)
	text(x(1),x(2)+0.01,x(3),'a_x','FontWeight','bold')
	text(y(1),y(2),y(3)+0.01,'a_y','FontWeight','bold')
	text(z(1)+0.01,z(2),z(3),'a_z','FontWeight','bold')
	
	%BPFZ lab frame sph tri
	pts2 = 1.2.*pts;
	t = num2cell(pts2,1);
	K = sphconvhulln(pts2,false);
	trisurf(K,t{:})
	title('lab frame')
	
	axis equal
	
	%BPFZ sample frame
	%normals
	nexttile
	t=num2cell(1.2*(R\nAlist.').',1);
	plot3(t{:},'r*')
	
	hold on
	
	t=num2cell(1.2*(R*nAlist.').',1);
	plot3(t{:},'b*')
	
	t=num2cell(1.2*nAlist,1);
	plot3(t{:},'c*')
	
	
	
	%---------true point--------
	t = num2cell(1.3*(R\nA.').',1); % (R\nA.').' == (R.'*nA.').'
	plot3(t{:},'ko')
	
	t = num2cell(1.3*(R.'\nA.').',1);
	plot3(t{:},'bo')
	
	t = num2cell(1.3*nA,1);
	plot3(t{:},'co')
	
	sphere(40);
	
	%BPFZ
	t = num2cell(1.2*(R\pts.').',1);
	K = sphconvhulln(pts,false);
	trisurf(K,t{:})
	title('sample frame')
	
	axis equal
	
	sgtitle(['geometry: ' geometry])
	
end

end



%%
%------------------------------CODE GRAVEYARD------------------------------
%{
%alternative to InCubicFZ.m
	%extract quaternions
	qlist = vertcat(symfive.q);
	qlist = disorientation(qlist,'cubic');

	%output vq points that intersect (and indices)
	nA_out = vertcat(five(inBPFZ).nA);


	inmisFZ = zeros(size(dlist,1),1);
	for j = 1:size(dlist,1)
		inmisFZ(i) = InCubicFZ(dlist);
	end

% 	q = five(i).q;
% 	nA = five(i).nA;

		% 		scl = 1.3;
		% 		x = A(:,1)*scl;
		% 		y = A(:,2)*scl;
		% 		z = A(:,3)*scl;
		%
		% 		w = [0,0,0].';
		% 		quiver3(w,w,w,x,y,z,0,'linewidth',1)
		% 		text(x(1),x(2)+0.01,x(3),'a_x','FontWeight','bold')
		% 		text(y(1),y(2),y(3)+0.01,'a_y','FontWeight','bold')
		% 		text(z(1)+0.01,z(2),z(3),'a_z','FontWeight','bold')

		% 		trisurf(K,t{:},'FaceColor','none')
		% 		xlim(xlims)
		% 		ylim(ylims)
		% 		zlim(zlims)
		
		% 		lims = [-1 1];
		% 		xlim(lims)
		% 		ylim(lims)
		% 		zlim(lims)

	
	% 	qtemp = disorientation(symfive(1).q,'cubic');
	% 	dtemp = q2rod(qtemp);
	
	% 	for j = 1:length(symfive)
	% 		symfive(i).q = qtemp;
	% 		symfive(i).d = dtemp;
	% 		symfive(i).geometry = geometry;

%}