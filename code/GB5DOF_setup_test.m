%GB5DOF_setup test
addpathdir({'cu2qu.m','q2rod.m','GBfive2oct.m'})
otmp = get_ocubo(1000);
o = get_octpairs(otmp);
five = GBoct2five(o);
propList = GB5DOF_setup(five);

figure
histogram(propList)
xlabel('GBE (J/m^2)')
ylabel('counts')