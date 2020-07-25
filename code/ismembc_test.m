%ismembc vs. ismember test

disp('long lists of numbers')
n=7*1e5; a=ceil(n*rand(n,1)); b=ceil(n*rand(n,1));

t1a = timeit(@() ismembc(a,b))
t1b = timeit(@() ismember(a,b))

disp('sorted long lists of numbers')
a = reshape(a,7,[]);
b = reshape(b,7,[]);

t2a = timeit(@() ismembc(a,b))
t2b = timeit(@() ismember(a,b))

disp('sparse array and set')
A = randi(floor(1e6/7),1e6,7);

%%
[locA,locB] = builtin('_ismemberhelper',sort(A),1:100);