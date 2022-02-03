function o = oflip(o)
arguments
   o(:,8) double = get_ocubo(2,'random',[],20)
end
% oflip  octonion flip, as in flip the rotation convention
o = [qinv(o(:,1:4)) qinv(o(:,5:8))];
end