function P = GaussElimModulo(A, b, q)
P = mod([A, b], q); 
[row_n, total_col] = size(P); 
m = size(A, 2); 

pivot_cols = zeros(row_n, 1); 
next_pivot_row = 1; 

for j = m : -1 : 1
    if next_pivot_row > row_n
        break; 
    end
    
    after_Pi = P(next_pivot_row:end, j);
    Index = find(after_Pi ~= 0);
    
    if isempty(Index)
        continue; 
    end
    
    actual_pivot_row = Index(1) + next_pivot_row - 1;
    P([next_pivot_row, actual_pivot_row], :) = P([actual_pivot_row, next_pivot_row], :);
    
    a = P(next_pivot_row, j);
    P(next_pivot_row, :) = mod(P(next_pivot_row, :) * minv(a, q), q);
    
    for k = 1:row_n
        if k ~= next_pivot_row && P(k, j) ~= 0
            factor = P(k, j);
            P(k, :) = mod(P(k, :) - factor * P(next_pivot_row, :), q);
        end
    end
    
    pivot_cols(next_pivot_row) = j;
    next_pivot_row = next_pivot_row + 1;
end

rank_r = next_pivot_row - 1;
if rank_r > 0
    valid_rows = P(1:rank_r, :);
    cols_found = pivot_cols(1:rank_r);
    
    [~, sort_idx] = sort(cols_found); 
    
    P_final = zeros(size(P));
    P_final(1:rank_r, :) = valid_rows(sort_idx, :);
    P = mod(P_final, q);
end