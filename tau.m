function simulated_t = tau(t0, target_population)
    % normal distribution
    t0 = 20;
    target_population = 300;

    rng("shuffle");
    % std = (max-min)/(2*Z)
    std = (30-1)/(2*1.96);

    %simulated_t = normrnd(t0,std,[1,target_population]);

    normal_dist = makedist('Normal', 'mu', t0, 'sigma', std);
    truncated_dist = truncate(normal_dist, 1, 30);
    simulated_t = round(random(truncated_dist, target_population, 1), 0);

    % skewed normal distribution
    % pearsrnd(mu,sigma,skew,kurt);
    % skew = 0.5;
    % kurt = 3;
    % [simulated_t,type] = pearsrnd(t0,std,skew,kurt,1, target_population);
    % disp(type);
    % simulated_t(simulated_t<1) = 1;

    figure;
    histogram(simulated_t, 'Normalization', 'count');
    xlabel('Tau (time horizon)');
    %ylabel('Population');
    %title('Simulated time horizon distribution among agents');
    xlim([0, 30]);
    ylim([1, 50]);
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
end