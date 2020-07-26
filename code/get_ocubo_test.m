%get_ocubo test
clear; close all;

test = 4;

seed = 10;
rng(seed)

addpathdir({'qinv.m'})

switch test
	case 1
		% generate a single octonion formed by two quaternions sampled randomly
		% from cubochoric space
		o = get_ocubo();
		
	case 2
		%generate 100 octonions randomly chosen from a list of quaternion
		%pairs generated using a cubochoric sidelength of 5 and uniform
		%sampling
		o = get_ocubo(100,'uniform',6);
		
	case 3
		%generate all combinations of quaternion pairs (i.e. octonions) using
		%5^3 uniformly sampled quaternions (7750 octonions), including
		%no-boundary octonions (2020-07-25, might specify  asoptional
		%name-value input argument later)
		o = get_ocubo([],'uniform',5);
		
	case 4
		% generate 100 random quaternion pairs (i.e. octonions)
		o = get_ocubo(100,'random');
		
	case 5
		% should produce an error since specifying both random cubochoric
		% sampling and a cubochoric grid sidelength is contradictory
		disp('this will produce an error')
		o = get_ocubo(100,'random',5);
end

disp(['# octonions generated: ' int2str(size(o,1))])

plotFZrodriguez_vtx();
five = GBoct2five(o);

ddisQ = true;
if ddisQ
	tmp = num2cell(q2rod(disorientation(vertcat(five.q),'cubic')),2);
	[five.d] = tmp{:};
end

% t = num2cell(q2rod(disorientation(vertcat(five.q),'cubic')),1);
% plot3(t{:},'*')
fig = figure;
fig.Position = [369.5000  416.5000  654.5000  292.5000];
tiledlayout(1,2)
plot5DOF(five,'ocubo')

plotmisFZtriQ = false;
if plotmisFZtriQ
	nexttile(1);
	d_all = vertcat(five.d);
	[r,c] = find(isinf(d_all));
	d_all(r,:)=[];
	K = convhulln(d_all);
	t = num2cell(d_all,1);
	trisurf(K,t{:})
end
