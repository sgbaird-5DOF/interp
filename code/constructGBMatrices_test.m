seed = 10;
rng(seed);
qA_Lab = randq();
qB_Lab = randq();
nA_Lab = normr(rand(3,1));
[gA_R,gB_R] = constructGBMatrices(qA_Lab,qB_Lab,nA_Lab,convention);