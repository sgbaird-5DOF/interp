% RANDQVALIDATE  validate randq.m
clear all; close all; clc
qtest = randq(1e6);
[w,t,p] = q2rot(qtest);

nw = 10;
nt = 10;
np = 20;
[edge_w1,edge_t1,edge_p1] = ndgrid(linspace(0,pi,nw),linspace(0,pi,nt),linspace(0,pi,np));
edge_w2 = edge_w1(2:end,2:end,2:end);
edge_t2 = edge_t1(2:end,2:end,2:end);
edge_p2 = edge_p1(2:end,2:end,2:end);
edge_w1 = edge_w1(1:end-1,1:end-1,1:end-1);
edge_t1 = edge_t1(1:end-1,1:end-1,1:end-1);
edge_p1 = edge_p1(1:end-1,1:end-1,1:end-1);

binvolume = 0.25.*(edge_p1-edge_p2).*(cos(edge_t1)-cos(edge_t2)).*(edge_w1-edge_w2-sin(edge_w1)+sin(edge_w2));

count = zeros(numel(edge_w1),1);
for i = 1:numel(count)
count(i) = sum(...
    edge_w1(i) <= w & w <= edge_w2(i) &...
    edge_t1(i) <= t & t <= edge_t2(i) &...
    edge_p1(i) <= p & p <= edge_p2(i));
end
count = count./binvolume(:);

figure;
bar(count)