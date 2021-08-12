function q = get_cubo(n,method,sidelength,printQ)
arguments
	n {mustBeNonNegIntegerOrEmpty} = 1
	method char {mustBeMember(method,{'random','uniform'})} = 'random'
	sidelength {mustBeNonNegIntegerOrEmpty} = double.empty
    printQ(1,1) logical = false
end
% GET_CUBO  get n quaternions from randomly or uniformly sampled cubochoric points
%--------------------------------------------------------------------------
% Inputs:
%		n - # of quaternions to output (re-calculated if using 'uniform' method
%		and sidelength is specified
%
%		method - sampling method, 'random', 'uniform'
%
%		sidelength - # of points along edge of cube used in cubochoric
%		sampling (automatically calculated if n is given and sidelength is
%		not specified)
%
% Outputs:
%		q - rows of quaternions
%
% Dependencies: cu2qu.m (De Graef CMU group, see GB Octonion code)
%
% References:
% [1] Singh, S., & De Graef, M. (2016). Orientation sampling for 
% dictionary-based diffraction pattern indexing methods. Modelling and
% Simulation in Materials Science and Engineering, 24(8).
% https://doi.org/10.1088/0965-0393/24/8/085013
%
% Author: Sterling Baird
%
% Date: 2020-07-25
%--------------------------------------------------------------------------
if strcmp(method,'uniform') && isempty(sidelength)
	sidelength = ceil(n^(1/3)); % auto-calculate sidelength
    n = sidelength^3;
elseif isempty(n)
	n = sidelength^3;
end

acube = 0.5*pi^(2/3);

%define grid
switch method
	case 'random'
		G = rand(n,3); %grid points
		
	case 'uniform'
		x1 = linspace(0,1,sidelength);
		[X,Y,Z] = ndgrid(x1,x1,x1);
		G = [X(:) Y(:) Z(:)]; %grid points
end

aa = acube*(2*G-1); %center grid about [0,0,0] and scale grid

% poolobj = gcp;
% addAttachedFiles(poolobj,{'cu2qu.m','cu2ho.m','ho2qu.m','ho2ax.m','ax2qu.m','GetPyramid.m'})

%convert to quaternion
q = zeros(n,4);
for i = 1:n
	q(i,:) = cu2qu(aa(i,:),printQ); %vectorization possible
end

end

%fortran code from https://github.com/EMsoft-org/EMsoft/GBs/EMBGO.f90
% acube = 0.5D0 * cPi**0.666666666
%aa = acube * (/ 2*rng_uniform(seed)-1.D0, 2.D0*rng_uniform(seed)-1.D0, 2.D0*rng_uniform(seed)-1.D0 /)
% qa = cu2qu(aa)

%----------------------CUSTOM VALIDATION FUNCTIONS-------------------------
function mustBeNonNegIntegerOrEmpty(arg)
errmsg = 'must be non-neg integer or empty []';
if ~isempty(arg)
	if floor(arg) ~= arg
		error(errmsg)
	elseif arg < 0
		error(errmsg)
	end
end
end
