e_s = 5.42181395e+01;
e_infty = 4;
tau = 2.81550596e-06;
alpha = 1 - 4.38316707e-01;

t_range = 10.^(-11:0.1:2);
z = t_range/tau;
% term1 = (e_s - e_infty) * (z .^ (alpha - 1)) / tau;
% term2 = mlf(alpha, alpha, -(z.^alpha), 15);
term1 = (e_s - e_infty) * (z .^ alpha);
term2 = mlf(alpha, alpha + 1, -(z .^ alpha), 20);

e = term1 .* term2;
semilogx(t_range, e);