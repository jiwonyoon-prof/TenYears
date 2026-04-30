function visualization_of_results(experimentType)

switch experimentType
    case 1
        visualize_ESS_and_RDD();
    case 2
        compare_TopK_with_varying_topK_and_qs();
    case 3
        compare_betamax_with_q();
end



function compare_betamax_with_q()
% Ensure the directory exists
if ~exist('./results', 'dir'), mkdir('./results'); end


fileName = './results/result_betamax_and_q.mat';

% check  whether the experiment has been done already.
if ~isfile(fileName)
    % if the file does not exist, generate the file by running actual code.
    calculate_betamax_with_varying_qs();
end
load('./results/result_betamax_and_q.mat')


figure(1); clf;

fSize = 22; 
set(gca, 'FontSize', fSize, 'LineWidth', 1.5);

mu = mean(beta_max, 2);
sigma = std(beta_max, 0, 2);

plot([setup.Primes(1), setup.Primes(end)], [setup.q, setup.q], 'LineWidth', 2.5, 'LineStyle', '-'); hold on;
plot([setup.Primes(1), setup.Primes(end)], [(2*setup.alpha+1)^setup.N, (2*setup.alpha+1)^setup.N], 'LineWidth', 2.5, 'LineStyle', '--');
errorbar(setup.Primes, mu, sigma, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'LineWidth', 2); 

grid on;

set(gca, 'YScale', 'log');
ytickformat('%.0f'); 

xlabel('$q_s$', 'FontSize', fSize+2, 'Interpreter', 'latex');
ylabel('Number of atomic operations', 'FontSize', fSize+2);


title_str = sprintf('$q=%d, k=%d, n=%d$', setup.q, setup.K, setup.N);
title(title_str, 'Interpreter', 'latex', 'FontSize', fSize+4);

lgd = legend(['$q=', num2str(setup.q), '$'], ...
             ['$(2\alpha+1)^n=', num2str((2*setup.alpha+1)^setup.N), '$'], ...
             '$\beta_{\max}$', ...
             'Interpreter', 'latex');
lgd.FontSize = fSize - 4;
lgd.Location = 'best';

save_path = ['./results/beta_max_and_qs(n=', num2str(setup.N),', q=', num2str(setup.q), '.png'];
print(gcf, save_path, '-dpng', '-r300');






function compare_TopK_with_varying_topK_and_qs()
setup= dosetup();
load(['./results/result_TopK(n=', num2str(setup.N), ', k=', num2str(setup.K), ').mat']);
figure(1); clf;
IDs = find(elapsedTimes(:, 4) == true);
notIDs = find(elapsedTimes(:, 4) == false);

plot3(elapsedTimes(IDs, 1), elapsedTimes(IDs, 2), elapsedTimes(IDs, 3), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
plot3(elapsedTimes(notIDs, 1), elapsedTimes(notIDs, 2), elapsedTimes(notIDs, 3), 'kx', 'MarkerSize', 8, 'LineWidth', 1.5);

grid on;

fSize = 15; 

xlabel('Top K', 'FontSize', fSize, 'FontWeight', 'bold');
ylabel('q_s', 'FontSize', fSize, 'FontWeight', 'bold');
zlabel('Elapsed Time (sec)', 'FontSize', fSize, 'FontWeight', 'bold');

set(gca, 'FontSize', 12);

lgd = legend('Succeeded to find the secret', 'Failed to find the secret');
set(lgd, 'FontSize', fSize, 'Location', 'northeast', 'Box', 'on'); 

fprintf('k=%d, n=%d, q=%d, \alpha=%f\n', setup.K, setup.N, setup.q, setup.alpha);
saveas(gcf,['./results/3d_comparison_for_topK(n=', num2str(setup.N),', q=', num2str(setup.q), '.png'])


figure(2); clf;
plot(elapsedTimes(IDs, 1), elapsedTimes(IDs, 3), 'ro', 'MarkerSize', 7); hold on;
plot(elapsedTimes(notIDs, 1), elapsedTimes(notIDs, 3), 'kx', 'MarkerSize', 7)
grid on;
xlabel('Top K');
ylabel('Elapased Time (sec)')
legend('Successed to find the secret!', 'Failed to find the secret');
fprintf('k=%d, n=%d, q=%d, \alpha=%f\n', setup.K, setup.N, setup.q, setup.alpha);
saveas(gcf,['./results/varingTopK(n=', num2str(setup.N),', q=', num2str(setup.q), '.png'])


figure(3); clf;
plot(elapsedTimes(IDs, 2), elapsedTimes(IDs, 3), 'ro', 'MarkerSize', 7); hold on;
plot(elapsedTimes(notIDs, 2), elapsedTimes(notIDs, 3), 'kx', 'MarkerSize', 7)
grid on;
xlabel('q_s');
ylabel('Elapased Time (sec)')
legend('Successed to find the secret!', 'Failed to find the secret');
fprintf('k=%d, n=%d, q=%d, \alpha=%f\n', setup.K, setup.N, setup.q, setup.alpha);
saveas(gcf,['./results/varing_gs(n=', num2str(setup.N),', q=', num2str(setup.q), '.png'])


figure(4); clf;
plot(elapsedTimes(IDs, 1), elapsedTimes(IDs, 2), 'ro', 'MarkerSize', 7); hold on;
plot(elapsedTimes(notIDs, 1), elapsedTimes(notIDs, 2), 'kx', 'MarkerSize', 7)
grid on;
xlabel('Top K');
ylabel('q_s');
legend('Successed to find the secret!', 'Failed to find the secret');
fprintf('k=%d, n=%d, q=%d, \alpha=%f\n', setup.K, setup.N, setup.q, setup.alpha);
saveas(gcf,['./results/varing_TopK_and_gs(n=', num2str(setup.N),', q=', num2str(setup.q), '.png'])



function visualize_ESS_and_RDD()
% Ensure the directory exists
if ~exist('./results', 'dir'), mkdir('./results'); end


fileName = './results/result_ESS_and_RDD.mat';

% check  whether the experiment has been done already.
if ~isfile(fileName)
    % if the file does not exist, generate the file by running actual code.
    compareEES_RDD();
end

% Load experimental data
load(fileName); % Expected: Ns, result4ESSRDD

num_Ns = length(Ns);
mu = zeros(num_Ns, 4);
sigma = zeros(num_Ns, 4);

for i = 1:num_Ns
    result = result4ESSRDD{i};
    % Ensure data is averaged correctly if multiple trials exist
    mu(i, :) = mean(result.elapasedTime, 1);
    sigma(i, :) = std(result.elapasedTime, 0, 1);
end

% Professional Color Palette
colors = [0, 0.4470, 0.7410;   % Blue
          0.8500, 0.3250, 0.0980;  % Red/Orange
          0.4660, 0.6740, 0.1880;  % Green (RDD Naive)
          0.4940, 0.1840, 0.5560]; % Purple (RDD Entropy)
          
labels = {'Conventional Error Enumeration (Naive)', 'Conventional Error Enumeration (Entropy)', ...
          'Proposed RDD (Naive)', 'Proposed RDD (Entropy)'};
lines = {'--', '--', '-', '-'};

fig1 = figure('Color', 'w', 'Position', [100, 100, 800, 600]); 
hold on;

% Prepare X-axis (ensure column vector for fill)
x_val = Ns(:); 

for i = 1:4
    y_val = mu(:, i);
    s_val = sigma(:, i);
    
    % 1. Plot Shaded Area (Standard Deviation)
    % Using 'HandleVisibility', 'off' so legends only show the mean lines
    fill_x = [x_val; flipud(x_val)];
    fill_y = [y_val - s_val; flipud(y_val + s_val)];
    fill(fill_x, fill_y, colors(i,:), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % 2. Plot Mean Line
    plot(x_val, y_val, lines{i}, 'LineWidth', 2.2, 'Color', colors(i,:), ...
        'DisplayName', labels{i});
        
    % 3. Add markers for clarity (optional)
    plot(x_val, y_val, 'o', 'MarkerSize', 4, 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'HandleVisibility', 'off');
end

% Formatting for Academic Publication
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlim([5, 10])
xlabel('Secret Dimension ($n$)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Execution Time (seconds)', 'FontSize', 16);
xticks(Ns);
%title('\textbf{Computational Efficiency: Error Enumeration vs. Proposed RDD}', 'Interpreter', 'latex', 'FontSize', 16);

% Move legend to the best location (usually northwest if complexity increases)
legend('Location', 'best', 'FontSize', 16);

% Save with high resolution
exportgraphics(fig1, './results/time_comparison_for_ESS_and_RDD.png', 'Resolution', 300);
fprintf('Plot saved to ./results/time_comparison_for_ESS_and_RDD.png\n');