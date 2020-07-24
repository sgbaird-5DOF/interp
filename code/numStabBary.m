function lambda = numStabBary(X,xA)
%{
------NOTES------

A numerically stable approach to computing barycentric coordinates by
converting to projective or homogeneous coordinates as an intermediate to
cartesian. I.e. cartesian -> projective -> barycentric. Uses an n-ary cross
product [1] computed with MATLAB's symbolic engine to solve a linear system
of equations (instead of using `\` (also called mldivide() or backsolve)
operator). Taking this approach should theoretically minimize division
operations that commonly cause numerical instability [2].

Max is 8D cartesian, add more symbolic vectors to U and I for higher
dimensions

------INPUT-------

X - row-wise vertices of simplex xA - test point that might be within
simplex

---INTERMEDIATE---

S - rows of vectors defining simplex, with X(1,:) as origin (i.e.
parametric representation)

------OUTPUT------

lambda - barycentric coordinates


----4D EXAMPLE----

X = [ 1 1 0 1; 1 1 1 0; 0 1 1 1; 1 0 1 1 ];

xA = [1 1 1 1];

lambda = numStabBary(X,xA);

----REFERENCES----

[1] https://en.wikipedia.org/wiki/Cross_product

[2] V. Skala, Robust Barycentric Coordinates Computation of the Closest
Point to a Hyperplane in E^n, Proc. 2013 Int. Conf. Applies Math. Comput.
Methods Eng. (2013) 239–244.

%}

d = length(xA);

x1 = X(1,:); %first vertex

S = zeros(size(X)-[1 0]); %one less row than X
srows = size(S,1);

% convert to parametric representation
for i = 2:size(X,1)
	S(i-1,:) = X(i,:)-x1;
end

% construct A matrix as in Ax=b
A = zeros(srows,srows);
for i = 1:srows
	for j = i:srows
		A(i,j) = dot(S(i,:),S(j,:)); %creates upper triangular matrix
	end
end
A = tril(A.') + triu(A,1);  % use upper triangle to make symmetric

%disp(['det(A) = ',num2str(det(A))]) %gives sense of stability (how thin the simplex is)

%construct b vector as in Ax=b
b = zeros(srows,1);
for i = 1:srows
	b(i) = dot(S(i,:),x1-xA);
end

method = 'extendedCross'; %'backslash', 'extendedCross'

switch method
	case 'extendedCross'
		%convert to extended cross product form (n-ary cross product)
		eta = zeros(srows,srows+1); %one more column to accomodate "b", one more row of 1's
		for i = 1:srows
			eta(i,:) = [A(i,:) b(i)];
		end
		
		%compute n-ary cross product (sigma)
		syms i1 i2 i3 i4 i5 i6 i7 i8 i9
		I = [i1 i2 i3 i4 i5 i6 i7 i8 i9]; %unit vectors
		eta = [eta; I(1:d)]; % Wikipedia suggests adding to last row so that direction of othorgonal vector is correct. Conflicts with reference in paper above
		etadet = det(eta); %can't just take determinant - top or bottom row needs to be i,j,k,l unit vectors, not 1,1,1,1
		
		if etadet == 0
			%disp('etadet == 0');
			lambda = repelem(-Inf,d);
			return
		end
		
		[C,T] = coeffs(etadet,I(1:d)); % k x 1 vector orthogonal to k-1 vectors
		Cvars = ismember(I(1:d),T);
		sigma(Cvars) = C;
		sigma(~Cvars) = 0;
		%compute parametric coordinates of point in Euclidean space
		
		U = double(sigma(1:end-1)/sigma(end));
	case 'backslash'
		U = A\b;
end

%compute barycentric coordinates
lambda = zeros(1,d);
lambda(1) = 1 - sum(U);

try
	lambda(2:end) = U;
catch
	lambda(2:end) = U;
end
%disp(lambda)
end %numStabBary


%{
Code Graveyard %eta(size(eta,1),:) = repelem(1,size(S,2)); %incorrect -
needs to be i,j,k,l unit vectors, not 1,1,1,1

etadet = det(eta); %can't just take determinant - top row needs to be i,j,k,l unit vectors, not 1,1,1,1

%}



