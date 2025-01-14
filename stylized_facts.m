function stylized_facts(log_returns)
    sample_size = 200;

    data = readtable('2022_2525_newmarket.csv');

    log_returns = data.LogReturn;
    %log_returns = log_returns(end-364:end);

    prices = data.Stock_Price;
    %prices = prices(end-364:end);
    figure;
    x = 1:length(prices);
    plot(x, prices);
    xlabel('Time Step');
    ylabel('Stock price changes');

    figure;
    x = 1:length(log_returns);
    plot(x, log_returns);
    xlabel('Time Step', 'FontSize', 12);
    ylabel('Log return of price', 'FontSize', 12);
    title('Time Series of Price Return', 'FontSize', 12);
    ylim([-0.03, 0.03]);
    %ylim([-0.1, 0.1]);

    set(findall(gcf,'-property','FontSize'),'FontSize',14)


    NormFitDemo(log_returns);
    hold on
    xlabel('Log return of price', 'FontSize', 12);
    title("Distribution of log price returns", 'FontSize', 12);
    xlim([-0.05, 0.05]);
    %xlim([-0.1, 0.1]);
    ylim([0, 120]);
    %ylim([0, 100]);
    hold off

    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    
    % auto_correlation
    lags_corr = 20;
    autocorr_val = autocorr(log_returns, 'NumLags', lags_corr);
    %autocorr = readtable('autocorr2022.csv');
    %autocorr = autocorr{:,:};
    %autocorr_val = mean(autocorr, 2);

    % calculate auto correlation
    figure;
    subplot(2,1,1);
    stem(0:lags_corr, autocorr_val(:,1), "filled");
    title('Auto-correlation of Log Return', 'FontSize', 12);
    xlabel('Lag', 'FontSize', 12);
    ylabel('Autocorrelation variable', 'FontSize', 12);
    hold on

    % 添加虚线边界
    yline(0); 
    alpha = 0.05; % critical interval
    %critical_value = norminv(1 - alpha / 2, 0, 1) / sqrt(sample_size);
    critical_value = 1.96 / sqrt(sample_size);
    yline(critical_value, '--');
    yline(-critical_value, '--');
    ylim([-0.5, 1]);
    
    %--------------------
    % absolute_correlation
    absolute_log_returns = abs(log_returns);
    autocorr_val_abs = autocorr(absolute_log_returns, 'NumLags', lags_corr);

    % calculate absolute correlation
    subplot(2,1,2);
    stem(0:lags_corr, autocorr_val_abs(:,1), "filled");
    title('Autocorrelation of absolute log return', 'FontSize', 12);
    xlabel('Lag', 'FontSize', 12);
    ylabel('Autocorrelation variable', 'FontSize', 12);

    hold on
    yline(critical_value, '--');
    yline(-critical_value, '--');
    ylim([-0.5, 1]);
    hold off

    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    
end