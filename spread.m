function spread = spread(bid_order_book, ask_order_book)

    if ~isempty(bid_order_book)
        bid_prices = [bid_order_book.price];
        bid_prices_levels = unique(bid_prices);
        best_bid = max(bid_prices_levels);
    else
        best_bid = NaN;
    end

    if ~isempty(ask_order_book)
        ask_prices = [ask_order_book.price];
        ask_prices_levels = unique(ask_prices);
        best_ask = min(ask_prices_levels);
    else
        best_ask = NaN;
    end

    spread = best_ask - best_bid;
    if isnan(spread)
        spread = 0.01;
    end
end