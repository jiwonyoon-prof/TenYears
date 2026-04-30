function [hat_hs, hat_b, D, hat_bs]= find_primitive_equation(setup, H, bs)
K = size(H, 1);
N = setup.N;

P = GaussElimModulo(H, bs, setup.q);
%P2 = modular_GaussElimination(H, bs, setup.q);


D = P(1:K-N, 1:N);
hat_hs = mod(P(1, 1:N+1)', setup.q);
hat_bs = mod(P(1:K-N, end), setup.q);
hat_b = hat_bs(1);

