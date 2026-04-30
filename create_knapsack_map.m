function [dp_val, dp_keep] = create_knapsack_map(setup, weights, counts, capacity, K)
    num_items = length(weights);
    dp_val = -inf(num_items + 1, capacity + 1, K);
    dp_keep = zeros(num_items + 1, capacity + 1, K, 2); 

    dp_val(1, 1, 1) = 0; % initial setting

    for i = 1:num_items
        w_i = weights(i);
        c_i = counts(i);

        for j = 0:capacity
            % All possible ways to choose $c$ units of the i-th item
            candidates_v = [];
            candidates_info = []; % [Selection Count, Parent $K$]

            for t = 0:c_i
                prev_cap = j - t * w_i;
                vt = (setup.alpha+10^(-10))^2/(2*setup.sigma^2)-(t-setup.alpha)^2/(2*setup.sigma^2);

                if prev_cap >= 0
                    prev_vals = squeeze(dp_val(i, prev_cap + 1, :));
                    valid_idx = find(prev_vals > -inf);
                    if ~isempty(valid_idx)
                        candidates_v = [candidates_v; prev_vals(valid_idx) + t * vt];
                        % For each value, store [Selection Count, Parent $K$]
                        for k_idx = 1:length(valid_idx)
                            candidates_info = [candidates_info; t, valid_idx(k_idx)];
                        end
                    end
                end
            end

            % Extracking the K best items
            if ~isempty(candidates_v)
                [sorted_v, sort_idx] = sort(candidates_v, 'descend');
                num_to_copy = min(length(sorted_v), K);
                dp_val(i+1, j+1, 1:num_to_copy) = sorted_v(1:num_to_copy);
                for k = 1:num_to_copy
                    dp_keep(i+1, j+1, k, :) = candidates_info(sort_idx(k), :);
                end
            end
        end
    end
end