% from Rodrigues vector to quaternions

function q = ro2qu(r,epsijk)
arguments
   r(:,4)
   epsijk(1,1) double = 1
end

npts = size(r,1);
q = zeros(npts,4);
for i = 1:npts
    rtmp = r(i,:);
    q(i,:) = ax2qu(ro2ax(rtmp),epsijk); %broken 2020-12-03
end
% set values very close to 0 as 0
thr = 1e-8;
q(abs(q) < thr) = 0;

% if (abs(q(1))-0)<thr
%     q(1)=0;
% elseif (abs(q(2))-0)<thr
%     q(2)=0;
% elseif (abs(q(3))-0)<thr
%     q(3)=0;
% elseif (abs(q(4))-0)<thr
%     q(4)=0;
% end