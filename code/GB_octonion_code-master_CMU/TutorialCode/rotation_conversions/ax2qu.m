% from axis-angle pair to quaternions
%vectorized by SGB 2020-08-15
function q = ax2qu(ax)

npts = size(ax,1);
q = zeros(npts,4);

thr = 1e-10;
ids = abs(ax(:,4)) < thr;

if any(ids)
	q(ids,:) = [ 1.0, 0.0, 0.0, 0.0 ];
end

if any(~ids)
	c = cos(ax(~ids,4)*0.5);
	s = sin(ax(~ids,4)*0.5);
	q(~ids,:) = [ c, ax(~ids,1).*s, ax(~ids,2).*s, ax(~ids,3).*s ];
end
% set values very close to 0 as 0

q(abs(q) < thr) = 0;

%--------------------------------CODE GRAVEYARD----------------------------
%{

thr = 1e-10;
if (abs(ax(4)-0.0)<thr)
   q = [ 1.0, 0.0, 0.0, 0.0 ];
else
   c = cos(ax(4)*0.5);
   s = sin(ax(4)*0.5);
   q = [ c, ax(1)*s, ax(2)*s, ax(3)*s ];
end

% set values very close to 0 as 0
if (abs(q(1))-0)<thr
    q(1)=0;
elseif (abs(q(2))-0)<thr
    q(2)=0;
elseif (abs(q(3))-0)<thr
    q(3)=0;
elseif (abs(q(4))-0)<thr
    q(4)=0;


%}