function [e_full, success] = enumeration_error_search(A, c, setup, e_1ton, gt_es)
[k, n]=size(A);
alpha = setup.alpha;
q = setup.q;

success = false;
e_full = []; s_rec = [];

total_candidates = size(e_1ton,1);

A1 = A(1:n, :);
c1 = c(1:n);
A2 = A(n+1:end, :);
c2 = c(n+1:end);

inv_A1 = inv_modulo(A1, q);
for i = 0 : total_candidates - 1
    e_part = e_1ton(i+1, :)';
    rhs = mod(c1 - e_part, q);
    try
        hat_s = mod(inv_A1 * rhs, q);
    catch
        continue; 
    end

    e_rest = allow_pm(mod(c2 - A2 * hat_s, q), q);
    
    if all(e_rest >= -alpha & e_rest <= alpha)
        e_full = [e_part; e_rest];
        success = true;
        return;
    end
end

