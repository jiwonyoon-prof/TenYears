function compareEES_RDD()
clc;
clear all;
close all;
% rand('seed', 3);
% rng('default');
rng(2);

setup.q = 48817;
setup.mu = 0;
setup.sigma = 1;
setup.alpha = 2*setup.sigma;

Ns = [5, 6, 7, 8, 9, 10];

generateErrors(setup, Ns);

nExperiments = 50;

for j=1:length(Ns)
    setup.N = Ns(j);
    setup.K = setup.N+2;
    load(['./results/full_errors', num2str(setup.N), '.mat']);

    fprintf('(n=%d, k=%d) for LWE is running.\n', setup.N, setup.K);
    for i=1:nExperiments
        fprintf('The %d th experiment is run.\n', i);

        %Generate LWE samples
        [cs, sk, A, es] = sample_from_LWE(setup);
        while length(find(es<-1*setup.alpha))>0||length(find(es>setup.alpha))>0
            [cs, sk, A, es] = sample_from_LWE(setup);
        end
  
        e_full = zeros(setup.K, 4);
  
        gt_es =es;

        myTime = cputime;
        [e_full(:, 1), success(1)] = enumeration_error_search(A, cs, setup, e_1ton_naive);
        elapasedTime(i, 1) =cputime-myTime;
        
        myTime = cputime;
        [e_full(:, 2), success(2)] = enumeration_error_search(A, cs, setup, e_1ton_entropy);
        elapasedTime(i, 2) =cputime-myTime;
        
        myTime = cputime;
        [e_full(:, 3), success(3)] = rdd_error_search(A, cs, setup, e_1ton_naive);
        elapasedTime(i, 3) =cputime-myTime;

        myTime = cputime;
        [e_full(:, 4), success(4)] = rdd_error_search(A, cs, setup, e_1ton_entropy);
        elapasedTime(i, 4) =cputime-myTime;
    
        hat_es{i} = e_full;
        if length(find(success==0))
            keyboard; % if there is no succeeded samples, then something is wrong!
        end
    end
    result4ESSRDD{j}.elapasedTime=elapasedTime;
    result4ESSRDD{j}.K = setup.K;
    result4ESSRDD{j}.N = setup.N;
    result4ESSRDD{j}.q = setup.q;
    result4ESSRDD{j}.alpha = setup.alpha;
    result4ESSRDD{j}.sigma = setup.sigma;
    result4ESSRDD{j}.nExperiments =nExperiments;
end

save('./results/result_ESS_and_RDD.mat', 'result4ESSRDD', 'Ns');




function generateErrors(setup, Ns)

% Setting error bounds: [-alpha, ..., 0, ..., alpha]
search_vals = -setup.alpha : setup.alpha;
range_size = length(search_vals);

for id4Ns=1:length(Ns)
    n = Ns(id4Ns);
    targetFile = ['./results/full_errors', num2str(n), '.mat'];
    if isfile(targetFile)
        continue;
    end

    total_candidates = range_size^n;
    
    e_1ton_naive = zeros(total_candidates, n);
    e_1ton_entropy = zeros(total_candidates, n);

    for i = 0 : total_candidates - 1
        temp_idx = i;
        for ith = n : -1 : 1
            val_idx = mod(temp_idx, range_size) + 1;
            % naive
            priority_vals = search_vals; % Naive: searching from -alpha
            e_1ton_naive(i+1, ith) = priority_vals(val_idx);
            % entropy
            [~, dist_idx] = sort(abs(search_vals));
            priority_vals = search_vals(dist_idx);
            e_1ton_entropy(i+1, ith) = priority_vals(val_idx);
            temp_idx = floor(temp_idx / range_size);
        end
    end
    save(targetFile, 'e_1ton_naive', 'e_1ton_entropy');
    fprintf('Full error vectors are made and saved for n=%d.\n', n);
end