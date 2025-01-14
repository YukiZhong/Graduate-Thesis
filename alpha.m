function simulated_a = alpha(a0, target_population)
    % normal distribution
    a0 = 0.3;
    target_population = 300;

    rng("shuffle");
    % std = (max-min)/(2*Z)
    std = 1/(2*1.96);
    
    %simulated_a = normrnd(a0,std,[1,target_population]);

    normal_dist = makedist('Normal', 'mu', a0, 'sigma', std);
    % truncate range limited to [0, 1]
    truncated_dist = truncate(normal_dist, 0, 1);
    simulated_a = random(truncated_dist, target_population, 1);

    % skewed normal distribution
    % pearsrnd(mu,sigma,skew,kurt);
    % skew = -0.5;
    % kurt = 3;
    % 
    % [simulated_a,type] = pearsrnd(a0,std,skew,kurt, 1, target_population);
    % disp(type);
    % simulated_a(simulated_a>1) = 1;
    
    figure;
    histogram(simulated_a, 'Normalization', 'count');
    %title('\rho Distribution of agents');
    xlabel('\rho');
    %ylabel('frequency');
    xlim([0, 1]);
    ylim([0, 100]);
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
end