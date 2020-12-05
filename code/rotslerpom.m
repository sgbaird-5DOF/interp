function r1 = rotslerpom(r0,q0,q1,dt,epsijk)
arguments
    r0(1,3) double = [1 0 0]
    q0(:,4) double = [1 0 0 0]
    q1(:,4) double = [0 0 0 1]
    dt(1,1) double = 0.1
    epsijk(1,1) double = 1
end
om = get_slerpom(q0,q1,dt,epsijk);
npts = length(om);
r1 = zeros(npts,3);
for i = 1:npts
%     r1(i,:) = (r0*om(:,:,i).').';
    r1(i,:) = (om(:,:,i)*r0.').';
%     r1(i,:) = (r0*(om{i}).').';
%     r1(i,:) = (om{i}\r0.').';
end
    

