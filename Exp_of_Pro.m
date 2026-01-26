d = 512;
alpha = 1.04;
q = 12289;
epsilon = 0.005;
beta = alpha - epsilon;
r = ((d - (d / 2) * beta^2) * q);
R = ((d / 2 * beta^2) * q);
lbound = r / (d / 24);
hbound = R / (d / 24);

b1 = sqrt(((d / 2 * alpha^2) * q) / (d / 24));
b2 = sqrt(((d - (d / 2) * alpha^2) * q) / (d / 24));

integrand = @(gamma) ((-marcumq(sqrt(gamma),b1,d) + marcumq(sqrt(gamma),b2,d)));

p = (1 / ((R / (d / 24)) - (r / (d / 24)))) * integral(integrand, lbound, hbound);