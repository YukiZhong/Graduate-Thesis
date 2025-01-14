function wealth_distribution()
    data = readtable('wealth_5_0.88_0.00012_newmarket.csv');
    data = table2array(data);
    disp(data);

    std_value = std(data);
    skewness_value = skewness(data);
    kurtosis_value = kurtosis(data);


    figure;
    bin_width = 500;
    %histogram(data, 'BinWidth', bin_width);  % set BinWidth
    histfit(data);
    xlim([-10000, 10000]);
    ylim([0, 50]);
    xlabel('Final wealth');
    ylabel('Frequency');
    title('Distribution of final wealth level');

    stats_text = sprintf('Std: %.2f\nSkewness: %.2f\nKurtosis: %.2f', ...
                     std_value, skewness_value, kurtosis_value);
    text(0.7 * max(data), 0.9 * max(histcounts(data, 'BinWidth', bin_width)), stats_text, ...
     'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'HorizontalAlignment', 'center');

    
    % sorted_data = data + abs(min(data));
    % sorted_data = sort(sorted_data, 'descend');
    % figure;
    % histogram(sorted_data, 'BinWidth', bin_width);

    %% log-log plot
    % figure;
    % loglog(1:length(sorted_data), sorted_data, 'bo-');
    % xlabel('Rank (log scale)');
    % ylabel('Final wealth (log scale)');
    % title('Log-Log Plot of Final Wealth');
end