function plot_spread_volume()
    spread = readtable('spreads_2022_2525_newmarket.csv');
    spread = table2array(spread);
    disp(length(spread));
    spread = spread(:, end-364:end);

    volume = readtable('volumes_2021_newmarket.csv');
    volume = table2array(volume);
    volume = volume(:, end-364:end);

    [num_runs, num_steps] = size(spread);

    time_steps = repmat(1:num_steps, num_runs, 1);
    time_steps = time_steps(:);
    spread_values = spread(:);
    volume_values = volume(:);

    % m = num_runs，n = num_steps
    %data = readtable('5_0.5_run5.csv');
    %avg_spread = data.avg_spread;
    %avg_volume = data.avg_volume;


    % figure;
    % boxchart(time_steps, spread_values, 'MarkerStyle', 'none');
    % ylim([-0.5, 4]);
    % 
    % xticks([1, num_steps]); 
    % xticklabels({'0', num_steps});
    % 
    % hold on;
    avg_spread = mean(spread, 1);
    % plot(1:num_steps, avg_spread, 'r-', 'LineWidth', 1);
    % hold off;
    % 
    % xlabel('Time Step');
    % ylabel('Spread');
    % title('Time Series Box Plot of Spread');
    % 
    % grid on;

    %%%%%%%%%%%%%%%%%
    avg_volume = mean(volume, 1); % mean for each time step
    
%%%%%%%%%%%%%%%
    
    % overall mean
    overall_spread_mean = mean(avg_spread);
    overall_volume_mean = mean(avg_volume);
    
    figure;
    bar(avg_spread);
    ylim([0, 2]);
    hold on;

    yline(overall_spread_mean, 'r-', 'LineWidth', 1.5);

    x_limits = xlim;
    y_limits = ylim;
    text(x_limits(1)+0.05*diff(x_limits), y_limits(2)-0.2*diff(y_limits), ...
            {['Mean = ', num2str(round(overall_spread_mean, 3))],['Std = ', num2str(round(std(avg_spread), 3))]}, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12);

    title('Time Series of Spread');
    xlabel('Time step');
    ylabel('Spread Between the Best Bid & Ask Prices');
    grid on;
    legend('Average of time series', 'Overall Average');


    %%%%%%%%%%%%%%%%%%%

    figure;
    bar(avg_volume); 
    ylim([0, 200000]);
    hold on;

    yline(overall_volume_mean, 'r-', 'LineWidth', 1.5);

    x_limits = xlim;
    y_limits = ylim;
    text(x_limits(1)+0.05*diff(x_limits), y_limits(2)-0.2*diff(y_limits), ...
            {['Mean = ', num2str(round(overall_volume_mean))],['Std = ', num2str(round(std(avg_volume)))]}, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12);

    title('Time Series of Volume');
    xlabel('Time step');
    ylabel('Transaction Volume');
    grid on;
    legend('Average of time series', 'Overall Average');


    disp(overall_volume_mean);
    disp(std(avg_volume));
    
    disp(overall_spread_mean);
    disp(std(avg_spread));

    
    figure;
    histogram(avg_spread, 'BinEdges', 0:0.01:2.5);
    xlim([0, 2.5]);
    ylim([0, 30]);

    title('Distribution of Bid–ask spread');
    xlabel('Spread Between the Best Bid & Ask Prices');
    ylabel('Frequency');

    
    NormFitDemo(avg_spread);
    hold on
    title("Distribution of Bid–ask spread");
    h.NumBins = 52;
    xlim([-0.1, 1.5]);
    ylim([0, 120]);
    hold off


end