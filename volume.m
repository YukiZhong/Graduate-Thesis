function [priceLevels, best_bid, best_ask] = volume(bid_order_book, ask_order_book)
    % disp("start");
    order_book = [ask_order_book, bid_order_book];
    
    best_bid = NaN;
    best_ask = NaN;

    if isempty(bid_order_book)
        bid_prices = [];
        bid_prices_levels = [];
    else
        bid_prices = [bid_order_book.price];
        bid_prices_levels = unique(bid_prices);
        best_bid = max(bid_prices_levels);
    end

    % disp("bid_prices");
    % disp(length(bid_order_book));
    
    if isempty(ask_order_book)
        ask_prices = [];
        ask_prices_levels = [];
    else
        ask_prices = [ask_order_book.price];
        ask_prices_levels = unique(ask_prices);
        best_ask = min(ask_prices_levels);
    end

    % disp("ask_prices");
    % disp(length(ask_order_book));

    priceLevels = struct('price', [], 'totalQuantity', []);
    all_prices = [order_book.price];
    all_levels = unique(all_prices);

    % disp("------------------");
    % disp(all_levels);
    % disp(best_bid);
    % disp(best_ask);

    for p = 1:length(all_levels)
        priceLevel = all_levels(p);
        idx = all_prices == priceLevel;
        totalQuantity = sum([order_book(idx).quantity]);
        
        priceLevels(p).price = priceLevel;
        priceLevels(p).totalQuantity = totalQuantity;
    end

    % for i = 1:length(priceLevels)
    %     disp(priceLevels(i));
    % end
end
    