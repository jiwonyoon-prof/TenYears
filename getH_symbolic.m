function [H, Hcs, newH, newb, est_es] =getH_symbolic(X, cs, q)

[k, n] = size(X);

XX = X'*X;


% P = modular_GaussElimination(XX, eye(n,n), q);
%P = GaussElimModulo(XX, eye(n,n), q);

inv_XX = new_modInv(XX, q);


H = mod(eye(k,k)-X*inv_XX*X', q);

HX = mod(H*X, q);

Hcs = mod(H*cs, q);

% Hes = mod(H*es, q);

% est_es = modular_GaussElimination(H, Hcs, q);

est_es = modGaussianElimination_symbolic(H, Hcs, q);

newH = est_es(1:k-n, 1:size(H, 2));
newb = est_es(1:k-n, end);

