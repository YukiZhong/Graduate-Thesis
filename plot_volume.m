function plot_volume(priceLevels, best_bid, best_ask)

    %prices = arrayfun(@(x) x.price, priceLevels);
    prices = [priceLevels.price];
    disp(prices);

    %quantities = arrayfun(@(x) x.totalQuantity, priceLevels);
    quantities = [priceLevels.totalQuantity];
    disp(quantities);
    
    maxPrice = max(prices);
    minPrice = min(prices);
    
    numBins = (maxPrice-minPrice)*100; 
    edges = linspace(minPrice, maxPrice, numBins + 4);

    [~, binIndex] = histc(prices, edges);
    binQuantities = zeros(1, numel(edges) - 1); 
    for i = 1:numel(edges) - 1
        binQuantities(i) = sum(quantities(binIndex == i));  
    end
    
    figure;
    hold on;
    
    %plot(edges, quantities, 'o', 'DisplayName', 'Price vs Quantity');
    bar(edges(1:end-1), binQuantities, 'histc');
    
    %x = edges;
    %y = ones(size(x)) * sum(quantities) / numel(prices);  % 假设数量沿着价格线均匀分布
    %plot(x, y, '--k', 'DisplayName', 'Price Range (Spread)');  % 绘制区间的 spread

    x_limits = xlim;
    y_limits = ylim; 
    bid_volume = sum(quantities(quantities > 0));
    ask_volume = sum(quantities(quantities < 0));


    if (~isnan(best_bid)) && (~isnan(best_ask))
        
        plot(best_bid, quantities(prices == best_bid), 'ro', 'MarkerSize', 10, 'DisplayName', 'Best Bid');
        plot(best_ask, quantities(prices == best_ask), 'go', 'MarkerSize', 10, 'DisplayName', 'Best Ask');
        
        % indicating spread
        spread = best_ask - best_bid;
        %text(mean([best_bid, best_ask]), max(quantities)*0.9, ['Spread: ', num2str(spread)], 'HorizontalAlignment', 'center', 'FontSize', 12);
        
        text(x_limits(1) + 0.05*diff(x_limits), y_limits(1) + 0.1*diff(y_limits), ...
            {['Spread: ', num2str(spread)],['Bid Volume: ', num2str(bid_volume)],['Ask Volume: ', num2str(ask_volume)]}, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12);
    else
        text(x_limits(1) + 0.05*diff(x_limits), y_limits(1) + 0.1*diff(y_limits), ...
            {['Bid Volume: ', num2str(bid_volume)],['Ask Volume: ', num2str(ask_volume)]}, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12);
    end
    
    xlabel('Price');
    ylabel('Total Quantity');
    title('Price vs Total Quantity');
    legend show;
    grid on;

    hold off;
end