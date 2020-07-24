function x = simplex_subdiv(v,n)
%---------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-06-30
%
% Description: Compute the subdivision cartesian coordinates of a simplex
%	using graded lexicographic compositions and equal subdivision of facet
%	edges. Basically just a wrapper function for John Burkardt's code.
%
% Input:
%
%		v		===	n x m coordinates of the vertices of the simplex where
%						m === dimension, i.e. rows of vertex coordinates
%
%		n		===	the number of subintervals (n == 1 does nothing, n == 2
%						does one subdivision, n == 3 does 2 subdivisions, etc.)
%
% Output:
%
%		x		===	the coordinates of one or more points (rows of points)
% 
% Note: adapted from John Burkardt code
% (http://people.math.sc.edu/Burkardt/m_src/simplex_grid/simplex_grid.html),
% which has various functions in the "HELPER FUNCTIONS" section.
%---------------------------------------------------------------------------

m = size(v,2); %dimension
ng = simplex_grid_size(m,n);
g = simplex_grid_index_all(m,n,ng);
x = simplex_grid_index_to_point(m,n,ng,g,v.');
x = x.';

end
%%
%{
%----------code graveyard---------
%x = x(1:d,:); %get rid of last row of zeros

%v = padarray(v,[0 1],'post');

%ng = simplex_grid_size(m,n);
%g = simplex_grid_index_all(m,n,ng);
%x = simplex_grid_index_to_point(m,n,ng,g,v);

%%
% count grid points inside a simplex
ng = 1;
for i = 1:d-1
	ng = (ng*(n+i))/i;
end
allcombQ = false;
%return simplex indices
if ~allcombQ
	grid = zeros(d+1,ng);
	
	g = zeros(d+1,1);
	g(d+1) = n;
	
	k = 1;
	grid(:,k) = g(:,1);
	
	for i = 1:ng
		g = comp_next_grlex(d+1,g);
		grid(:,i) = g(:,1);
	end
	
else
	vec = 0:n-1;
	mat = repelem({vec},d);
	grid = allcomb(mat{:});
	grid = sortrows([sum(grid,2) grid]);
end

nQ = true;
if ~nQ
	take = grid(:,1) == d-1;
else
	take = 1:ng;
end

%grid = grid(take,2:end).';

%%
%return points corresponding to simplex indices
x = (v.'*grid/n).';

% Dependencies:
%	No longer needed (2020-06-30)
%		allcomb.m, MATLAB FEX,
%		https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin
%}

%%
%-------------------------HELPER FUNCTIONS---------------------------------
function xc = comp_next_grlex(kc,xc)

%*****************************************************************************80
%
%% COMP_NEXT_GRLEX returns the next composition in grlex order.
%
%  Discussion:
%
%    Example:
%
%    D = 3
%
%    #   XC(1) XC(2) XC(3)  Degree
%      +------------------------
%    1 |  0     0     0        0
%      |
%    2 |  0     0     1        1
%    3 |  0     1     0        1
%    4 |  1     0     0        1
%      |
%    5 |  0     0     2        2
%    6 |  0     1     1        2
%    7 |  0     2     0        2
%    8 |  1     0     1        2
%    9 |  1     1     0        2
%   10 |  2     0     0        2
%      |
%   11 |  0     0     3        3
%   12 |  0     1     2        3
%   13 |  0     2     1        3
%   14 |  0     3     0        3
%   15 |  1     0     2        3
%   16 |  1     1     1        3
%   17 |  1     2     0        3
%   18 |  2     0     1        3
%   19 |  2     1     0        3
%   20 |  3     0     0        3
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 December 2013
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer KC, the number of parts of the composition.
%
%    Input, integer XC(KC), the current composition.
%    The lowest order composition is XC = [ 0, 0, ..., 0, 0 ].
%
%    Output, integer XC(KC), the next composition.
%

%
%  Ensure that 1 <= KC.
%
  if ( kc < 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'COMP_NEXT_GRLEX - Fatal error!' );
    fprintf ( 1, '  KC < 1\n' );
    error ( 'COMP_NEXT_GRLEX - Fatal error!' );
  end
%
%  Ensure that 0 <= XC(I).
%
  for i = 1 : kc
    if ( xc(i) < 0 )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'COMP_NEXT_GRLEX - Fatal error!' );
      fprintf ( 1, '  XC(I) < 0\n' );
      error ( 'COMP_NEXT_GRLEX - Fatal error!' );
    end
  end
%
%  Find I, the index of the rightmost nonzero entry of XC.
%
  i = 0;
  for j = kc : -1 : 1
    if ( 0 < xc(j) )
      i = j;
      break
    end
  end    
%
%  set T = XC(I)
%  set XC(I) to zero,
%  increase XC(I-1) by 1,
%  increment XC(KC) by T-1.
%
  if ( i == 0 )
    xc(kc) = 1;
    return
  elseif ( i == 1 )
    t = xc(1) + 1;
    im1 = kc;
  elseif ( 1 < i )
    t = xc(i);
    im1 = i - 1;
  end

  xc(i) = 0;
  xc(im1) = xc(im1) + 1;
  xc(kc) = xc(kc) + t - 1;

  return
end


%%
%{
%-----------------------------CODE GRAVEYARD-------------------------------
% while (k < ng)
% 	g = comp_next_grlex(m+1,g);
% 	k = k + 1;
% 	grid(:,k) = g(:,1);
% end

%}

%%
%--------------------------HELPER FUNCTION GRAVEYARD-----------------------
%%
function ng = simplex_grid_size ( m, n )

%*****************************************************************************80
%
%% SIMPLEX_GRID_SIZE counts the grid points inside a simplex.
%
%  Discussion:
%
%    The size of a grid with parameters M, N is C(M+N,N).
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    30 July 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of subintervals.
%
%    Output, integer NG, the number of grid points.
%
  ng = 1;
  for i = 1 : m
    ng = ( ng * ( n + i ) ) / i;
  end

  return
end

%%
function grid = simplex_grid_index_all ( m, n, ng )

%*****************************************************************************80
%
%% SIMPLEX_GRID_INDEX_ALL returns all the simplex grid indices.
%
%  Discussion:
%
%    The number of grid indices can be determined by calling 
%      ng = simplex_grid_size ( m, n );
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    30 July 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of subintervals.
%
%    Input, integer NG, the number of values in the grid.
%
%    Output, integer GRID(M+1,NG), the current, and then the next, grid index.
%
  grid = zeros ( m + 1, ng );

  g = zeros ( m + 1, 1 );
  g(m+1,1) = n;

  k = 1;
  grid(1:m+1,k) = g(1:m+1,1);

  while ( k < ng )
    g = comp_next_grlex ( m + 1, g );
    k = k + 1;
    grid(1:m+1,k) = g(1:m+1,1);
  end

  return
end




function x = simplex_grid_index_to_point( m, n, ng, g, v )

%*****************************************************************************80
%
%% SIMPLEX_GRID_INDEX_TO_POINT returns  points corresponding to simplex indices.
%
%  Discussion:
%
%    The M-dimensional simplex is defined by M+1 vertices.
%
%    Given a regular grid that uses N subintervals along the edge between
%    each pair of vertices, a simplex grid index G is a set of M+1 values
%    each between 0 and N, and summing to N. 
%
%    This function determines the coordinates X of the point corresponding
%    to the index G.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    31 July 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of subintervals.
%
%    Input, integer NG, the number of grid indices to be converted.
%
%    Input, integer G(M+1,NG), the grid indices of 1 or more points.
%
%    Input, real V(M,M+1), the coordinates of the vertices of the simplex.
%
%    Output, real X(M,NG), the coordinates of one or more points.
%
  x = v(1:m,1:m+1) * g(1:m+1,1:ng) / n;
	return
end