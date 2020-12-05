function om = get_slerpom(q0,q1,dt,epsijk)
arguments
    q0(:,4) double = [1 0 0 0]
    q1(:,4) double = [0 0 0 1]
    dt(1,1) double = 0.1
    epsijk(1,1) double = 1
end
vq = get_slerp(q0,q1,dt);
om = qu2om(vq,epsijk);

end

%% CODE GRAVEYARD
%{
% npts = size(vq,1);
% om = cell(npts,1);
% for i = 1:npts
%     vqtmp = vq(i,:);
%     om{i} = qu2om(vqtmp,epsijk);
% end
%}