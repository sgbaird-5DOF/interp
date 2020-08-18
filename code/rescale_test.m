%rescale test (test of built-in function)

clear; close all

%generate data
seed = 10;
rng(10);
data = [0;abs(lognrnd(10,1,100000,1))];

%histograms
fig = figure;
fig.Position = [355.8 139 380 420];
nexttile
histogram(data,'LineStyle','none')
hold on
histogram(-data+max(data),'LineStyle','none')
legend('data','-data+max(data)','Location','best')
ylabel('counts')
xlabel('property')

nexttile
histogram(rescale(data),'LineStyle','none')
hold on
histogram(rescale(-data),'LineStyle','none')

legend('rescale(data)','rescale(-data)','Location','best')
ylabel('counts')
xlabel('property')