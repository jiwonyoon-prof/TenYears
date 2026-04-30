function main()

clc;
clear all;
close all;
rand('seed', 3);
rng('default');
% rng(2);

setup= dosetup();



%Generate LWE samples
[cs, sk, A, es] = sample_from_LWE(setup);
while length(find(es<-1*setup.alpha))>0||length(find(es>setup.alpha))>0
    [cs, sk, A, es] = sample_from_LWE(setup);
end


%[as_validation, cs_validation, es_validation] = getValidationData(setup, sk)

% A, cs, sk, es



sk, es
gt_es = es; % only for simulation for comparing with the ground truth....

debugging = false;
elapsedTimes = zeros(length(setup.multipleTopK)*length(setup.Primes), 4);
for j=1:length(setup.multipleTopK)
    setup.TopK = setup.multipleTopK(j);
    fprintf('\n\n Top K (=%d) is selected.\n', setup.TopK);
    for i=1:length(setup.Primes)
        qs = setup.Primes(i);
        fprintf('The used q_s is %d\n', qs);
        timePoint = cputime;
        [hat_es, isSolved, ws_mod_qs] = Solve_LWE_Search(setup, cs, A, qs);
        elapsedTimes((j-1)*length(setup.Primes)+i, :) = [setup.TopK, qs, cputime-timePoint, isSolved];

        check_qs(setup.q, qs, ws_mod_qs)
    end
end
elapsedTimes


if ~isfolder('./results')
    mkdir('./results');
end
save(['./results/result_TopK(n=', num2str(setup.N), ', k=', num2str(setup.K), ').mat'], 'elapsedTimes', 'setup');


visualization_of_results(1);
visualization_of_results(2);
visualization_of_results(3);



function check_qs(q, qs, ws_mod_qs)
if length(ws_mod_qs)~= length(unique(ws_mod_qs))
    disp('==> Note that redundant value is detected in ws_mod_qs! This makes myriads of unwanted local opima.')
end

if length(find(ws_mod_qs==0))>0
    disp('==> Note that current qs makes zero value in ws^*. This makes myriads of unwanted local opima')
end




function [hat_es, isSolved, ws_mod_qs] = Solve_LWE_Search(setup, cs, A, qs)

K = setup.K;
N = setup.N;

hat_es = zeros(K, 1)+999;

S = [];
isSolved = false;

% Build H only with K-1 data
[H, Hcs] =getH(A, cs, setup.q);    

% Build bs
bs = Hcs + setup.alpha*H*ones(size(A, 1), 1);  % bs = H*ts  where bs = H+alpha*H*1



% Build Aumented matrix
% Find primitive equation, ws'*ts(1:N+1) = b
[ws, hat_b, D, hat_bs]= find_primitive_equation(setup, H, bs);



% changing real to mod q_s
ws_mod_qs = mod(ws, qs); %  % ({\bf w})_{q_s} = {\bf w} mod q_s


beta_max = ceil(2*setup.alpha*sum(ws_mod_qs)); 
r_max = ceil( (ceil(2*setup.alpha)*sum(ws)-hat_b)/setup.q );
ws_star = ws_mod_qs;

[dp_val, dp_keep] = create_knapsack_map(setup, ws_star, 2*setup.alpha*ones(1, N+1), beta_max, setup.TopK);

full_ks =[];

hat_ts = zeros(K, 1);
for r=0:r_max
    b_1_mod_qs_plus_rq_with_fixed_r = mod(hat_b+r*setup.q, qs);
    f_max_with_r = ceil( (ceil(2*setup.alpha) *sum(ws_mod_qs) - b_1_mod_qs_plus_rq_with_fixed_r)/qs );
    for f=0:f_max_with_r

        beta = b_1_mod_qs_plus_rq_with_fixed_r + f*qs;
        if beta>beta_max
            % disp('beta cannot be larger than beta_{max}.');
            continue;
        end
        Trajs = [];
        % [Trajs, ~] = getTopKSet(M, P, R, ws_star, N+1, beta);
        for k=1:setup.TopK
            final_val = dp_val(end, beta+1, k);
            if final_val<0
                continue;
            end
            full_ks = [full_ks, k];
            items = trace_knapsack_items(dp_keep, ws_star, beta, k);
            Trajs = [Trajs, items'];
        end
        if size(Trajs, 2)==0
            continue
        end

        nTrajs = size(Trajs, 2);

        hat_ts_rears=mod(repmat(hat_bs, 1, nTrajs) - D*Trajs(1:N, :), setup.q);
        temp_ts_rears = allow_pm(hat_ts_rears, setup.q);
        myLogic = sum((temp_ts_rears>=0) & (temp_ts_rears<= (2*setup.alpha)), 1);
        ID4es = find(myLogic==2);
        if length(ID4es) > 0
            hat_ts(1:N) = Trajs(1:N, ID4es(1)); 
            hat_ts(N+1:K) = hat_ts_rears(:, ID4es(1));
            hat_es = hat_ts - setup.alpha*ones(K, 1);
            isSolved = true;
            return;
        end
    end
end
