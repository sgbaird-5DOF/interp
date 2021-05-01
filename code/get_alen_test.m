%get_alen_test
seed = 10;
rng(seed);
a = normr(rand(1,3)); %generate random point on 2-sphere
b = normr(rand(2,3)); %generate two random points on 2-sphere
alen = get_alen(a,b)