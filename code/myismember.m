function ids = myismember(a,b)

prec = 6;
r = @(a) round(a,prec);
tol = 1e-3;
ids = ismembertol(r(a),r(b),tol,'ByRows',true);

end