%get_ocubo test
clear; close all;

test = 3;

seed = 10;
rng(seed)

addpathdir({'qinv.m'})

switch test
	case 1
		% generate a single octonion formed by two quaternions sampled randomly
		% from cubochoric space
		o = get_ocubo();
		
	case 2
		% generate 5 random octonions
		o = get_ocubo(5,'random');
		
	case 3
		%generate ~100 uniform octonions
		o = get_ocubo(100,'uniform');
		disp(['# octonions generated: ' int2str(size(o,1))])
		
	case 4
		%generate all combinations of quaternion pairs (i.e. octonions) using 5^3
		%uniformly sampled quaternions (7750 octonions), including no-boundary
		%octonions (2020-07-25, might specify optional name-value input argument
		%later)
		o = get_ocubo([],'uniform',5);
		
	case 5
		% generate 100 random quaternion pairs (i.e. octonions) sampled from 5^3
		% uniformly sampled quaternions
		o = get_ocubo(100,'random',5);
		
	case 6
		% generate 100 quaternion pairs (i.e. octonions) from randomly sampled
		% quaternions
		o = get_ocubo(100,'random');
end

plotFZrodriguez_vtx();
five = GBoct2five(o);

t = num2cell(q2rod(disorientation(vertcat(five.q),'cubic')),1);
plot3(t{:},'*')
