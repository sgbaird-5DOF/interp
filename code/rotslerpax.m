function [x1,y1,z1] = rotslerpax(x,y,z,q0,q1,dt,epsijk)
arguments
    x(1,3) double = [1 0 0]
    y(1,3) double = [0 1 0]
    z(1,3) double = [0 0 1]
    q0(:,4) double = [1 0 0 0]
    q1(:,4) double = [0 0 0 1]
    dt(1,1) double = 0.1
    epsijk(1,1) double = 1
end

% x1 = rotslerpom(x,q0,q1,dt,epsijk);
% y1 = rotslerpom(y,q0,q1,dt,epsijk);
% z1 = rotslerpom(z,q0,q1,dt,epsijk);

x1 = rotslerp(x,q0,q1,dt,epsijk);
y1 = rotslerp(y,q0,q1,dt,epsijk);
z1 = rotslerp(z,q0,q1,dt,epsijk);

end