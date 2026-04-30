function P = modGaussianElimination_symbolic(A, b, q)

[n, m] = size(A);
q_sym = sym(q);

aug = mod([sym(A), sym(b)], q_sym);

target_row = n;

for j = m:-1:1
    if target_row < 1, break; end
    
    found = false;
    for r = target_row:-1:1
        if mod(aug(r, j), q_sym) ~= 0
            actual_row = r;
            found = true;
            break;
        end
    end
    
    if found
        if actual_row ~= target_row
            aug([target_row, actual_row], :) = aug([actual_row, target_row], :);
        end
        
        pivot = mod(aug(target_row, j), q_sym);
        [~, invP, ~] = gcd(pivot, q_sym);
        invP = mod(invP, q_sym);
        aug(target_row, :) = mod(aug(target_row, :) * invP, q_sym);
        
        for i = 1:n
            if i ~= target_row
                factor = mod(aug(i, j), q_sym);
                aug(i, :) = mod(aug(i, :) - factor * aug(target_row, :), q_sym);
            end
        end
        
        target_row = target_row - 1;
    end
end

isZeroRow = all(aug == 0, 2);
aug(isZeroRow, :) = [];

P = aug;