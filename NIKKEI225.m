function NIKKEI225()

    data = readtable('HistoricalData_SP500.csv');
    %data = sortrows(data(1008:end, :), 'Date'); %2020
    %data = sortrows(data(756:1007, :), 'Date'); %2021

    data = readtable('NIKKEI225.csv');
    data.DATE = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');
    disp(data);
    %data = data(~isnan(data.Close), :);

    data = data(~isnan(data.NIKKEI225), :);

    % log return time series
    figure;
    plot(data.DATE, data.NIKKEI225);
    plot(data.DATE, data.Close);
    title('Daily Stock Average');
    xlabel('Date');
    ylabel('Price');
    grid on;

    logReturns = {};
    whole_prices = [];

    bootstrapSamples = 1000;
    fat_tails = zeros(bootstrapSamples, 1);
    autocorr_log_returns = [zeros(bootstrapSamples, 20)];
    autocorr_abs_log_returns = [zeros(bootstrapSamples, 20)];

    logReturns_all = log(data.NIKKEI225(2:end) ./ data.NIKKEI225(1:end-1));
    dataTable0 = table(data.NIKKEI225(2:end), logReturns_all, 'VariableNames', {'Stock_Price', 'LogReturn'});

    %logReturns_all = log(data.Close(2:end) ./ data.Close(1:end-1));
    %dataTable0 = table(data.Close(2:end), logReturns_all, 'VariableNames', {'Stock_Price', 'LogReturn'});
    writetable(dataTable0, 'Nikkei2020.csv');

    % for currentyear = 2019:2019
    % 
    %     currentYearData = data(year(data.DATE) == currentyear, :);
    %     disp(currentYearData);  
    % 
    %     %prices = currentYearData.NIKKEI225;
    %     whole_prices = [whole_prices; prices];
    % 
    %     currentYearReturn = log(prices(2:end) ./ prices(1:end-1));
    % 
    %     n = length(currentYearReturn);
    %     bootstrapReturns = datasample(currentYearReturn, n * bootstrapSamples, 'Replace', true);  
    %     reshapedReturns = reshape(bootstrapReturns, n, bootstrapSamples);
    % 
    %     %disp(reshapedReturns);
    % 
    %     logReturns = [logReturns; currentYearReturn];
    % 
    % end
end