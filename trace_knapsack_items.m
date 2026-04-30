function item_counts = trace_knapsack_items(dp_keep, weights, target_cap, k)
num_items = length(weights);
item_counts = zeros(1, num_items);

curr_cap = target_cap;
curr_k = k;

for i = num_items + 1 : -1 : 2
    info = squeeze(dp_keep(i, curr_cap + 1, curr_k, :));
    num_selected = info(1);
    parent_k = info(2);
    
    item_counts(i-1) = num_selected;
    
    curr_cap = curr_cap - num_selected * weights(i-1);
    curr_k = parent_k;
    
    if curr_k == 0, break; end
end
