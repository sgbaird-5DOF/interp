%structvertcat_test

S1 = struct('var1',1,'var3',{'hello'})
S2 = struct('var1',2,'var2',[4 5 6],'var4',struct())

Sout = structvertcat(S1,S2)
Sout(1)
Sout(2)