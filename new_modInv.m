function A_inv = new_modInv(A, q)

A_sym = sym(A);

d = det(A_sym);
d_mod = mod(d, q);

if d_mod == 0
    error('No inverse matrix!');
end

[~, n1, ~] = gcd(sym(d_mod), q);
d_inv = mod(n1, q);

adjA = adjoint(A_sym);
A_inv = mod(adjA * d_inv, q);

