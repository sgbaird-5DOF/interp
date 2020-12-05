function ax = vecplot(scl,r,w,mkr,rlbl)
arguments
   scl(1,1) double = 1
   r(1,3) double = [1,0,0]
   w(1,3) double = [0,0,0]
   mkr char = 'k'
   rlbl = 'r';
end
% VECPLOT  plot three quivers in the x-hat, y-hat, and z-hat directions
r = scl*r;
wtmp = n2c(w);
t = n2c(r);
hold on
ax(1) = quiver3(wtmp{:},t{:},0,mkr,'linewidth',1,'Autoscale','off');
ax(2) = text(r(1),r(2)+0.05,r(3),rlbl,'FontWeight','bold');

end

%% CODE GRAVEYARD
%{
% x = x.'*scl;
% y = y.'*scl;
% z = z.'*scl;
%}