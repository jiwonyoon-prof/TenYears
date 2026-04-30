function [cs, vec_s, vec_a, e] = sample_from_LWE(setup)
q = setup.q;
n = setup.N;
k = setup.K;

vec_s = mod(round(rand(n, 1)*setup.q^2), q); % secret
vec_a = mod(round(rand(k, n)*setup.q^2), q); % public
%e = round(normrnd(0, (alpha*q), k, 1)); % e~\chi : that is, e is generated from a discrete normal distribution.
e = mod(round(normrnd(setup.mu, setup.sigma, k, 1)), q); % e~\chi : that is, e is generated from a discrete normal distribution.
cs = mod(vec_a*vec_s + e, q);


e = allow_pm(e, q);