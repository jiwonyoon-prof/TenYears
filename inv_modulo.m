function A_inv = inv_modulo(A, p)
sz = size(A, 1);
Aug = mod([A, eye(sz)], p);
for i = 1:sz
    [~, inv_val, ~] = gcd(double(Aug(i,i)), p);
    Aug(i,:) = mod(Aug(i,:) * inv_val, p);
    for j = 1:sz
        if i ~= j
            factor = Aug(j,i);
            Aug(j,:) = mod(Aug(j,:) - factor * Aug(i,:), p);
        end
    end
end
A_inv = Aug(:, sz+1:end);
