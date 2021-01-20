function A = sigfn(x,scl,xshift)
arguments
    x
    scl(1,1) double = 30
    xshift(1,1) double = 1.1
end
% SIGFN  sigmoid function coefficients
A = 1./(1+exp(-scl*(x-xshift)));
end