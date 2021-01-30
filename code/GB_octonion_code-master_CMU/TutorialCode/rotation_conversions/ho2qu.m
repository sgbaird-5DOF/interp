% from homochoric to quaternions

function q = ho2qu(h,printQ)
arguments
    h
    printQ(1,1) logical = false
end
if printQ
    disp('ho2ax')
end
ax = ho2ax(h);
if printQ
    disp('ax2qu')
end
q = ax2qu(ax);

% set values very close to 0 as 0
thr = 1e-8;
if (abs(q(1))-0)<thr
    q(1)=0;
elseif (abs(q(2))-0)<thr
    q(2)=0;
elseif (abs(q(3))-0)<thr
    q(3)=0;
elseif (abs(q(4))-0)<thr
    q(4)=0;
end