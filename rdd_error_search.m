function [e_full, success] = rdd_error_search(A, cs, setup, e_1ton)
[k,n]=size(A);
q= setup.q;
alpha = setup.alpha;

success = false;
e_full = [];

total_candidates = size(e_1ton, 1);

[H, bs] =getH(A, cs, setup.q);
[~, ~, DH, b_prime]= find_primitive_equation(setup, H, bs);


for i = 0 : total_candidates - 1
    e_part = e_1ton(i+1, :)';%zeros(n, 1);
    e_rest = mod(b_prime - DH * e_part, q);
    e_rest = allow_pm(e_rest, q);
    if all(e_rest >= -alpha & e_rest <= alpha)
        e_candidate = [e_part; e_rest];
        is_valid = true;
        e_full = e_candidate;
        success = true;
        return;
    else
        e_candidate = [];
        is_valid = false;
    end
end

