function calculate_betamax_with_varying_qs()

setup.N = 8; % this is for the dimension of secret
setup.K = setup.N+2; % k dimension for A ==> k*n
% setup.multipleTopK = [5, 10, 50, 100, 200, 300, 400, 500];
setup.q = 4782971;  % Use for large q.
setup.mu = 0;%(setup.q-1)/2;
setup.sigma = 1;
setup.small_val= 10^(-10);
setup.alpha = 2*setup.sigma;

FullPrimes = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997];
ID = sort(unique(round(rand(1, 20)*length(FullPrimes))));
setup.Primes = FullPrimes(ID);
Ind = find(setup.Primes>setup.N*2);
setup.Primes = setup.Primes(Ind(1):end);


q = setup.q;

fprintf('q=%d, k=%d, n=%d', q, setup.K, setup.N);

N = 20; % the number of experiments

for i=1:length(setup.Primes)
    qs = setup.Primes(i);
    fprintf('The used q_s is %d\n', qs);
    for j=1:N
        [cs, sk, A, es] = sample_from_LWE(setup);
        while length(find(es<-1*setup.alpha))>0||length(find(es>setup.alpha))>0
            [cs, sk, A, es] = sample_from_LWE(setup);
        end
        [H, Hcs] =getH_symbolic(A, cs, q);    
        
        % Build bs
        bs = Hcs + setup.alpha*H*ones(size(A, 1), 1);  % bs = H*ts  where bs = H+alpha*H*1
        
        % Build Aumented matrix
        % Find primitive equation, ws'*ts(1:N+1) = b
        [ws, hat_b, D, hat_bs]= find_primitive_equation_symbolic(setup, H, bs);
        
        % changing real to mod q_s
        ws_mod_qs = mod(ws, qs); %  % ({\bf w})_{q_s} = {\bf w} mod q_s
        
        beta_max(i, j) = ceil(2*setup.alpha*sum(ws_mod_qs)); 
    end
end
save('./results/result_betamax_and_q.mat', 'beta_max', 'setup');
