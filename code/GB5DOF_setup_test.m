%GB5DOF_setup test
otmp = get_ocubo(1000);
o = get_octpairs(otmp);
[qm, nA] = oct2five(o);
propList = GB5DOF_setup([], qm, nA);
propList2 = GB5DOF_setup(o(:,1:4), o(:,5:8));

figure
histogram(propList)
xlabel('GBE (J/m^2)')
ylabel('counts')