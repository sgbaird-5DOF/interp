%VECPAIR2RMAT_TEST  vecpair2rmat gives proper rotations that take v1 to v2
% Covers antipodal pairs, for which vecpair2rmat used to return -eye(3)
% (an improper rotation that om2qu turns into a complex quaternion), as
% well as identical and generic pairs.
%
% There are no %% sections on purpose: run-tests runs each section as a
% separate test in its own workspace.

rng(1)
v = normr(randn(4,3));
w = normr(randn(4,3));

% antipodal (incl. the octonion convention's BP normal), identical, generic
V1 = [0 0 1; 1 0 0; v; v; v];
V2 = [0 0 -1; -1 0 0; -v; v; w];

tol = 1e-9;
for i = 1:size(V1,1)
    v1 = V1(i,:);
    v2 = V2(i,:);

    % active (default)
    R = vecpair2rmat(v1,v2,1);
    assert(abs(det(R)-1) < tol,'det(R) = %g for pair %d',det(R),i)
    assert(norm(R.'*R-eye(3)) < tol,'R is not orthogonal for pair %d',i)
    assert(norm(R*v1.'-v2.') < tol,'R*v1 ~= v2 for pair %d',i)
    assert(isreal(om2qu(R,1)),'om2qu gives a complex quaternion for pair %d',i)

    % passive
    R = vecpair2rmat(v1,v2,-1);
    assert(abs(det(R)-1) < tol,'det(R) = %g for pair %d (epsijk = -1)',det(R),i)
    assert(norm(R.'*v1.'-v2.') < tol,'R.''*v1 ~= v2 for pair %d (epsijk = -1)',i)
end

% a GB whose BP normal is antipodal to the octonion convention's [0 0 1]
o = five2oct([1 0 0 0],[0 0 -1]);
assert(isreal(o) && abs(norm(o)-sqrt(2)) < tol,'five2oct fails for nA = [0 0 -1]')
