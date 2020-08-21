function q = qu2ax(qq)
%--------------------------------------------------------------------------
% Date: 2020-08-21
%
% Description: from quaternions to axis-angle pair
% 
% Inputs:
%  qq - rows of quaternions
%
% Outputs:
%  q - rows of axis-angle pairs
%
% Usage:
%  q = qu2ax(qq);
%
% Notes:
%  Vectorized by SGB (2020-08-21)
%--------------------------------------------------------------------------
thr = 1e-8;
omega = 2*acos(qq(:,1));

ids = omega < thr;
if any(ids)
    q(ids,:) = [0, 0, 1, 0];
end

ids2 = ~ids & (qq(:,1) < thr);
if any(ids2)
    q(ids2,:) = [qq(2), qq(3), qq(4), pi];
end

if any(~ids2)
    s = (qq(~ids2,1)./abs(qq(~ids2,1)))./rssq(qq(~ids2,2:4),2);
    q(~ids2,:) = [qq(~ids2,2:4).*[s s s], omega(~ids2,:)];
end

% set values very close to 0 as 0
q(q < thr) = 0;

end



%---------------------------CODE GRAVEYARD---------------------------------
%{
% s =  (qq(~ids2,1)/abs(qq(~ids2,1)))/sqrt(qq(~ids2,2).^2+qq(~ids2,3).^2+qq(~ids2,4).^2);

%------Original CMU group function------
function q = qu2ax(qq)

thr = 1e-8;
omega = 2.0 * acos(qq(1));

if ((omega-0.0)<thr)
  q = [0.0, 0.0, 1.0, 0.0];
else
    
  if ((qq(1)-0.0)<thr)
      q = [qq(2), qq(3), qq(4), pi];
  else
    s =  (qq(1)/abs(qq(1)))/sqrt(qq(2)^2+qq(3)^2+qq(4)^2);
    q = [ qq(2)*s, qq(3)*s, qq(4)*s, omega];
  end
 
end

% set values very close to 0 as 0
if (abs(q(1))-0)<thr
    q(1)=0;
elseif (abs(q(2))-0)<thr
    q(2)=0;
elseif (abs(q(3))-0)<thr
    q(3)=0;
end
%}