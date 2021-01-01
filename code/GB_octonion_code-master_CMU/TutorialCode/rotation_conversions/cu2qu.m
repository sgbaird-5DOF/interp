% from cubochoric to quaternions

function q = cu2qu(c)
disp('cu2ho')
h = cu2ho(c);
disp('ho2qu')
q = ho2qu(h);

% set values very close to 0 as 0
thr = 1e-10;
if (abs(q(1))-0)<thr
    q(1)=0;
elseif (abs(q(2))-0)<thr
    q(2)=0;
elseif (abs(q(3))-0)<thr
    q(3)=0;
elseif (abs(q(4))-0)<thr
    q(4)=0;
end