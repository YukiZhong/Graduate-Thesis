% use “Command” + “/” to comment and “Command” + “Option” + “/” to uncomment.

clear;

% Model parameters
%num_steps = 3950;
num_steps = 2125;
initial_stock_price = 100;
market_order_ratio = [0, 0.1, 0.3, 0.5, 0.7, 0.9]; % <- modify
num_agents = 200;

% Define the types of investors and their settings
investor_types = {'Random', 'CARA', 'IARA'};
initial_wealth_IARA = 30;

% Preallocate arrays to store results
num_runs = 5;

% Call the main function to run the simulation
agent_based_model(num_runs, investor_types, market_order_ratio, num_steps, initial_stock_price, num_agents, initial_wealth_IARA);


% Main function definition
function agent_based_model(num_runs, investor_types, market_order_ratio, num_steps, initial_stock_price, num_agents, initial_wealth_IARA)

    % -------- plot fundamental values ------------
    fundamental_values2 = zeros(num_steps+1, 1);
    mu_f = 0.01;
    sigma_f = 0.05;
    dt_f = 1/1000;

    rng default;

    % Generate the predicted price of the stock
    sigma_N = 0.5;

    fundamental_values2(1, 1) = initial_stock_price;
    for t = 1:num_steps
        fundamental_values2(t+1, 1) = fundamental_values2(t, 1) * exp(dt_f*(mu_f-0.5*sigma_f^2)+sigma_f*sqrt(dt_f)*normrnd(0,1,1));
    end
    
    
    % plot fundamental values from skip+
    skip = 300;
    unit = 5;
    figure;
    x=(1:(num_steps-skip)/unit)';
    plot(x, fundamental_values2(skip+2:unit:end));
    title('Fundamental value (information)');
    xlabel('Time Steps');
    ylabel('Fundamental values of stock price');
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    %----------------------------------------------

    num_people_each_option = NaN;
    t2_list = NaN;
    sample_size = NaN;

    xi = NaN;
    f = NaN;

    for i = 2:length(investor_types)-1
        for j = 3:length(market_order_ratio)-3
    
            % initialize
            stock_prices_runs = zeros(num_runs, num_steps+1);
            price_volume = cell(num_runs, 1);

            volumes_runs = zeros(num_runs, (num_steps-skip)/unit+1);
            spreads_runs = zeros(num_runs, (num_steps-skip)/unit+1);

            daily = zeros(num_runs, (num_steps-skip)/unit+1);
            log_return_runs = zeros(num_runs, (num_steps-skip)/unit);
            wealths = zeros(num_runs, num_agents);

            % ----------------------------------------------
            
            % can change to other years
            startDate = datetime(2022, 1, 1);
            dates = startDate + days(0:(num_steps - skip)/unit-1);
            
            dataTable1 = table();
            dataTable2 = table();
            for run = 1:num_runs
                fprintf("run %d \n", run);
                
                % Initialise order_book
                bid_order_book = [];
                ask_order_book = [];
                quantity_unit = 0;

                % calculate for each agent; and for all agents
                %executed_orders = [];
    
                current_price = initial_stock_price; % Initialise current_price
                stock_prices_runs(run, 1) = initial_stock_price;
                
                [tau_distribution, num_people_each_option, t2_list, sample_size] = Time_pref(num_agents, num_people_each_option, t2_list, sample_size);
                [alpha_distribution, xi, f] = estimate_rho(num_agents, xi, f);
                
                %tau_distribution = tau(5, num_agents);
                %alpha_distribution = alpha(0.88, num_agents);


                %initialize agents
                investors=[];
                for k = 1:num_agents
                    
                    mu = 0;
            
                    if i == 3
                        initial_wealth = initial_wealth_IARA; % Initialise initial_wealth
                        wealth = 30 * ones(1, num_steps);
                    else
                        initial_wealth = 0;
                        wealth = zeros(1, num_steps);
                    end
                    
                    tau_0 = tau_distribution(k)*unit;
                    alpha_0 = alpha_distribution(k);
                    investors = [investors, me_agents(k, mu, initial_wealth, wealth, investor_types(i), num_steps, initial_stock_price, tau_0, alpha_0)];
                end
   
                % for each time step
                % to get the distribution of quantity
                quantity_t = [];

                for t = 1:num_steps
                    temporary_order_book = [];
                    
                    % for each agent, set details
                    for k = 1:num_agents
                        investors(k) = investors(k).cal_time(t);

                        % Set weights of fundamental value, price history and noise to random values (0, 1)
                        n_F = rand; % Following uniform distribution U(0,1)
                        n_P = rand; % Following uniform distribution U(0,1)
                        n_N = rand; % Following uniform distribution U(0,1)
                        investors(k) = investors(k).calculate_tau(n_F, n_P, n_N);

                        % Calculate the average price over the past memory_length time steps
                        investors(k) = investors(k).calculate_p_bar(stock_prices_runs(run, :), t);

                        % Calculate the basic value of the stock using geometric Brownian motion
                        p_F = fundamental_values2(t, 1);

                        % Generate the noisy part
                        epsilon = randn * sigma_N; % Following normal distribution N(0, sigma_N^2)
                        
                        investors(k) = investors(k).generate_predicted_price(p_F, epsilon, current_price);

        
                        % at the beginning, conduct call auction to determine opening price
                        if t == 1
        
                            % Decide the parameters for the order (if applicable)
                            investors(k) = investors(k).decide_order_type(market_order_ratio(j), true);

                            investors(k) = investors(k).quote_price(bid_order_book, ask_order_book, current_price);

                            % calculate the new position for the investor
                            investors(k) = investors(k).calculate_position(t, stock_prices_runs(run, :));
        
                            % make orders only when quantity is not 0
                            if investors(k).quantity_t ~= 0
                                new_order = investors(k).create_orders(t);
                                
                                % Add the new orders directly to the order book
                                if strcmp(new_order.bs_type, 'buy')
                                    bid_order_book = add_to_book(bid_order_book, new_order);
                                    
                                else
                                    ask_order_book = add_to_book(ask_order_book, new_order);
                                end
                            end
                        end
                    end

                    % execute
                    % call auction at the beginning

                    if t == 1
                        [bid_order_book, ask_order_book, opening_price, investors] = call_auction(bid_order_book, ask_order_book, current_price, investors);

                        current_price = opening_price;
                        disp("opening_price");
                        disp(opening_price);

                        % plot the volume
                        %[priceLevels, best_bid, best_ask] = volume(bid_order_book, ask_order_book);
                        %[priceLevels, best_bid, best_ask] = volume(bid_order_book, ask_order_book);
                        %plot_volume(priceLevels, best_bid, best_ask);

                    end

                    for k = 1:num_agents

                        % Calculate the basic value of the stock using geometric Brownian motion
                        p_F = fundamental_values2(t+1, 1);

                        % Generate the predicted price of the stock
                        % Following normal distribution N(0, sigma_N^2)
                        epsilon = randn * sigma_N;
                        investors(k) = investors(k).generate_predicted_price(p_F, epsilon, current_price);
 
                        % Add the new order to the temporary order book
                        % Decide the parameters for the order (if applicable)
                        investors(k) = investors(k).decide_order_type(market_order_ratio(j), false);

                        investors(k) = investors(k).quote_price(bid_order_book, ask_order_book, current_price);

                        % calculate the new position for the investor
                        investors(k) = investors(k).calculate_position(t, stock_prices_runs(run, :));
                        current_quantity = investors(k).quantity_t;

                        quantity_t = [quantity_t, current_quantity];


                        % make orders only when quantity is not 0
                        if investors(k).quantity_t ~= 0
                            new_order = investors(k).create_orders(t);
                            temporary_order_book = [temporary_order_book, new_order];
                        end
                    end

                    % for the rest of the times
                    % Sort the order book by time
                    if ~isempty(temporary_order_book)
                        temporary_order_book = sort_orders_by_time(temporary_order_book);
                
                        % Execute the orders in the order book
                        [bid_order_book, ask_order_book, new_price, investors, quantity] = stock_exchange(temporary_order_book, bid_order_book, ask_order_book, current_price, investors);
                        quantity_unit = quantity_unit + quantity;

                        current_price = new_price;
                    end

                    %%% end of a trading day
                    stock_prices_runs(run, t+1) = current_price;
                    
                    % average return of all agents at each time step t
                    %sum_utility = 0;

                    for k = 1:num_agents
                        % Calculate the wealth and utiltiy after transactions in this time step
                        investors(k) = investors(k).calculate_wealth_t(t, current_price);
   
                        %sum_utility = sum_utility + investors(k).utility_ts(t);
                    end

                    %[priceLevels, best_bid, best_ask] = volume(bid_order_book, ask_order_book);
                    %plot_volume(priceLevels, best_bid, best_ask);

                    % disp("=======show=======");
                    % if isempty(bid_order_book)
                    %     for p = 1:5
                    %         disp(ask_order_book(end-p));
                    %     end
                    % else
                    %     for p = 1:length(bid_order_book)
                    %         disp(bid_order_book(p));
                    %         disp(ask_order_book(end-p+1));
                    %     end
                    % end

                    
                    % if t>10
                    %     return
                    % end

                    
                    if(rem(t,unit) == 0)
                        fprintf("day %d \n", t/unit);
                        disp(current_price);
                        if (t >= skip)
                            daily(run, (t-skip)/unit+1) = current_price;

                            %[priceLevels, best_bid, best_ask] = volume(bid_order_book, ask_order_book);
                            %plot_volume(priceLevels, best_bid, best_ask);
                            %order_book_volumes{run, (t-skip)/unit+1} = priceLevels;
                            
                            current_spread = spread(bid_order_book, ask_order_book);
                            spreads_runs(run, (t-skip)/unit+1) = current_spread;
                            if isnan(current_spread)
                                error("spread is nan");
                            end
                            volumes_runs(run, (t-skip)/unit+1) = quantity_unit;
                        end
                        
                        % automatic diminishing of unmatched orders
                        % cancel = unit*30;
                        % if t > cancel
                        %     bid_order_book = bid_order_book((t - [bid_order_book.time]) <= cancel);
                        %     ask_order_book = ask_order_book((t - [ask_order_book.time]) <= cancel);
                        % end

                        quantity_unit = 0;
                    end
                end

                [priceLevels, best_bid, best_ask] = volume(bid_order_book, ask_order_book);
                plot_volume(priceLevels, best_bid, best_ask);

                figure;
                histogram(quantity_t);                

                %disp("executed_orders");
                %disp(length(executed_orders));
               
                for k = 1:num_agents
                    wealths(run, k) = investors(k).cal_final_profit();
                end

                log_returns = [];
                for a = 2:length(daily(run, :))
                    next = log(daily(run, a)/daily(run, a-1));
                    log_returns = [log_returns,next];
                end

                log_return_runs(run, :) = log_returns;

                [priceLevels, best_bid, best_ask] = volume(bid_order_book, ask_order_book);
                price_volume{run} = priceLevels;

            end
            
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % after all simulation runs, plot the figures

            % !!!!average stock price
            avg_stock_prices = sum(daily, 1)/num_runs;
            
            dataTable = table(dates', avg_stock_prices(1,2:end)', 'VariableNames', {'Date', 'Stock_Price'});
            return_from_avg_price = [];
            for a = 2:length(avg_stock_prices)
                next = log(avg_stock_prices(a)/avg_stock_prices(a-1));
                return_from_avg_price = [return_from_avg_price,next];
            end
            dataTable.return_from_avg_price = return_from_avg_price';

            %average log return
            avg_returns = sum(log_return_runs, 1)/num_runs;
            dataTable.LogReturn = avg_returns';
            disp(dataTable);

            avg_spread_runs = round(sum(spreads_runs, 1)/num_runs, 2);
            dataTable.avg_spread = avg_spread_runs(1,2:end)';

            avg_volume_runs = sum(volumes_runs, 1)/num_runs;
            dataTable.avg_volume = avg_volume_runs(1,2:end)';

            avg_wealths = sum(wealths, 1)/num_runs;
            %disp("sum of wealths for each run");
            %disp(avg_wealths);
            
            writetable(dataTable, '2022_2525_newmarket.csv');
            writematrix(avg_wealths, 'wealth_2022_2525_newmarket.csv');

            writematrix(volumes_runs, 'volumes_2022_2525_newmarket.csv');
            writematrix(spreads_runs, 'spreads_2022_2525_newmarket.csv');


            %!!! average order book volume
            all_levels = struct('price', [], 'totalQuantity', []);

            for run = 1:num_runs
                % extract priceLevels for current time step
                pricelevel = price_volume{run};
                if ~isempty(pricelevel)
                    for p = 1:length(pricelevel)
                        current_level = pricelevel(p);

                        if isempty(all_levels)
                            all_levels(1) = current_level;
                        else
                            price_pos = find([all_levels.price] == current_level.price, 1);
                            if isempty(price_pos) || price_pos <= 0
                                all_levels = [all_levels; current_level];
                            else
                                all_levels(price_pos).totalQuantity = all_levels(price_pos).totalQuantity + current_level.totalQuantity;
                            end
                        end
                    end

                end
            end


            % create bin

            [~, order] = sort([all_levels.price]);
            sorted_levels = all_levels(order);
            disp(sorted_levels);
            plot_volume(sorted_levels, NaN, NaN);

        end

        % for each agent type
        %figure;
        %plot(1:length(market_order_ratio), std_price_changes);
        %title("agent type " + investor_types(i));
        %xlabel('market_order_ratio');
        %ylabel('standard deviation of price change');
    end

    % print overall average total profit
    %for i = 1:length(investor_types)
        %for j = 1:length(market_order_ratio)
            %fprintf('average_total_profit_per_group(%d, %d) = %f\n', i, j, total_avg_profits(i, j));
        %end
    %end
end



% modified function
function [new_bid_order_book, new_ask_order_book, new_price, investors, quantity] = stock_exchange(temporary_order_book, bid_order_book, ask_order_book, current_price, investors)
    % set default value of new_price/new_position to current_price/current_position
    new_ask_order_book = ask_order_book;
    new_bid_order_book = bid_order_book;
    orders = temporary_order_book;
    quantity = 0;

    %disp("length of temporary order book");
    %disp(length(orders));

    % execute orders
    for i = 1:length(orders)
        order = orders(i);
        %[new_ask_order_book, new_bid_order_book, new_executed_orders, new_price, investors] = match_orders(new_ask_order_book, new_bid_order_book, order, current_price, investors);
        [new_ask_order_book, new_bid_order_book, new_price, investors, new_quantity] = match_orders(new_ask_order_book, new_bid_order_book, order, current_price, investors);
        quantity = quantity + new_quantity;
        
        current_price = new_price;
        %executed_orders = [executed_orders, new_executed_orders];
    end
end
    


% modified to include execution of a single temporary order
function [new_ask_order_book, new_bid_order_book, new_price, investors, quantity] = match_orders(ask_order_book, bid_order_book, order, current_price, investors)
    % initialise
    %new_executed_orders = [];
    new_ask_order_book = ask_order_book;
    new_bid_order_book = bid_order_book;

    n = 0;  % market order executed times
    exchange_price = current_price;
    new_price = exchange_price;
    quantity = 0;


    if(strcmp(order.bs_type, 'buy'))
        for i = length(new_ask_order_book):-1:1
            
            % define prices first
            % consider market orders
            if strcmp(order.type, 'market')
                % reach execute limit
                if (n >= 5)
                    order.type = 'limit';
                    order.price = current_price;
                    new_bid_order_book = add_to_book(new_bid_order_book, order);
                    return
                end

                n = n+1;
                
                exchanged_temp_order = order;
                opposite_order = new_ask_order_book(i);
                exchanged_oppo_order = opposite_order;
                
                if strcmp(opposite_order.type, 'market')
                    % ! 假设每个time step内也更新价格，使得market order达成时使用新价格
                    %exchange_price = current_price;
                    exchanged_oppo_order.price = exchange_price;
                    exchanged_temp_order.price = exchange_price;
                else
                    exchange_price = opposite_order.price;
                    exchanged_temp_order.price = exchange_price;
                end

            % if opposite is market
            elseif strcmp(new_ask_order_book(i).type, 'market')
                exchanged_temp_order = order;
                opposite_order = new_ask_order_book(i);
                exchanged_oppo_order = opposite_order;

                exchange_price = order.price;
                exchanged_oppo_order.price = exchange_price;

            else
                % if both are limit order, check if there is a match
                if is_matching(order, new_ask_order_book(i))
    
                    exchanged_temp_order = order;
                    opposite_order = new_ask_order_book(i);
                    exchanged_oppo_order = opposite_order;
                    
                    % exchange price is first come price
                    % if order.time < opposite_order.time
                    %     exchange_price = order.price;
                    %     exchanged_oppo_order.price = exchange_price;
                    % else
                    %     exchange_price = opposite_order.price;
                    %     exchanged_temp_order.price = exchange_price;
                    % end
                    exchange_price = opposite_order.price;
                    exchanged_temp_order.price = exchange_price;

                else
                    % no matching, means there is a gap between prices
                    break
                end
            end
            
            new_price = exchange_price;

            
            % consider quantity if there are matching orders
            %if ~isnan(exchange_price)

            % if temporary order has smaller quantity
            if abs(order.quantity) <= abs(opposite_order.quantity)
                exchanged_oppo_order.quantity = -order.quantity;
                quantity = quantity + abs(order.quantity);
                
                if abs(order.quantity) == abs(opposite_order.quantity)
                    new_ask_order_book(i) = [];
                   
                else
                    % !
                    new_ask_order_book(i).quantity = opposite_order.quantity + order.quantity;
                end

                % @
                investors(exchanged_oppo_order.agent) = investors(exchanged_oppo_order.agent).record_update(exchanged_oppo_order);
                investors(exchanged_temp_order.agent) = investors(exchanged_temp_order.agent).record_update(exchanged_temp_order);
                %new_executed_orders = [new_executed_orders, exchanged_temp_order, exchanged_oppo_order];
                
                % current temporary order executed
                return

            elseif abs(order.quantity) > abs(opposite_order.quantity)
                % !
                exchanged_temp_order.quantity = -opposite_order.quantity;
                order.quantity = order.quantity - exchanged_temp_order.quantity;
                quantity = quantity + abs(opposite_order.quantity);

                new_ask_order_book(i) = [];
                %new_executed_orders = [new_executed_orders, exchanged_temp_order, exchanged_oppo_order];
                % @
                investors(exchanged_oppo_order.agent) = investors(exchanged_oppo_order.agent).record_update(exchanged_oppo_order);
                investors(exchanged_temp_order.agent) = investors(exchanged_temp_order.agent).record_update(exchanged_temp_order);
            end

        end

        % resume the loop
        new_bid_order_book = add_to_book(new_bid_order_book, order);
        return
    
    elseif (strcmp(order.bs_type, 'sell'))

        for i = length(new_bid_order_book):-1:1
            % define prices first
            % consider market orders
            if strcmp(order.type, 'market')
                % reach execute limit
                if (n >= 5)
                    order.type = 'limit';
                    order.price = current_price;
                    new_ask_order_book = add_to_book(new_ask_order_book, order);
                    return
                end

                n = n+1;
                
                exchanged_temp_order = order;
                opposite_order = new_bid_order_book(i);
                exchanged_oppo_order = opposite_order;
                
                if strcmp(opposite_order.type, 'market')
                    %exchange_price = current_price;
                    exchanged_oppo_order.price = exchange_price;
                    exchanged_temp_order.price = exchange_price;
                else
                    exchange_price = opposite_order.price;
                    exchanged_temp_order.price = exchange_price;
                end

            % if opposite is market
            elseif strcmp(new_bid_order_book(i).type, 'market')
                exchanged_temp_order = order;
                opposite_order = new_bid_order_book(i);
                exchanged_oppo_order = opposite_order;

                exchange_price = order.price;
                exchanged_oppo_order.price = exchange_price;

            else
                % if both are limit order, check if there is a match
                if is_matching(order, new_bid_order_book(i))
    
                    exchanged_temp_order = order;
                    opposite_order = new_bid_order_book(i);
                    exchanged_oppo_order = opposite_order;
                    
                    % exchange price is first come price

                    exchange_price = opposite_order.price;
                    exchanged_temp_order.price = exchange_price;

                else
                    % no matching, means there is a gap between prices
                    break
                end
            end
            
            new_price = exchange_price;


            % if temporary order has smaller quantity
            if abs(order.quantity) <= abs(opposite_order.quantity)
                exchanged_oppo_order.quantity = -order.quantity;
                quantity = quantity + abs(order.quantity);

                if abs(order.quantity) == abs(opposite_order.quantity)
                    new_bid_order_book(i) = [];
                else
                    new_bid_order_book(i).quantity = opposite_order.quantity + order.quantity;
                end

                %new_executed_orders = [new_executed_orders, exchanged_temp_order, exchanged_oppo_order];
                investors(exchanged_oppo_order.agent) = investors(exchanged_oppo_order.agent).record_update(exchanged_oppo_order);
                investors(exchanged_temp_order.agent) = investors(exchanged_temp_order.agent).record_update(exchanged_temp_order);
                
                % current temporary order executed
                return

            elseif abs(order.quantity) > abs(opposite_order.quantity)
                exchanged_temp_order.quantity = -opposite_order.quantity;
                order.quantity = order.quantity - exchanged_temp_order.quantity;
                quantity = quantity + abs(opposite_order.quantity);
                new_bid_order_book(i) = [];

                %new_executed_orders = [new_executed_orders, exchanged_temp_order, exchanged_oppo_order];
                investors(exchanged_oppo_order.agent) = investors(exchanged_oppo_order.agent).record_update(exchanged_oppo_order);
                investors(exchanged_temp_order.agent) = investors(exchanged_temp_order.agent).record_update(exchanged_temp_order); 
            end
        end
        % resume the loop
        new_ask_order_book = add_to_book(new_ask_order_book, order);
    end
end


% see if there is at least one match
function is_match = is_matching(order1, order2)
    if (strcmp(order1.bs_type, 'buy') && order1.price >= order2.price)
        is_match = true;
    elseif (strcmp(order1.bs_type, 'sell') && order1.price <= order2.price)
        is_match = true;
    else
        is_match = false;
    end
end


% add uncomplete temporary orders into order book
function new_order_book = add_to_book(new_order_book, order)

    % !!
    if strcmp(order.type, 'market')
        % find matching index
        index = find(strcmp({new_order_book.type},order.type));
    else
        index = find_price_index(order.price, new_order_book);
    end

    if ~isempty(index)
        % if order with a same price level already exist, add the new order
        % by time
        % use descending time
        insert_pos = index(find([new_order_book(index).response_time] < order.response_time, 1, 'first'));
        
        if isempty(insert_pos)
            insert_pos = index(end) + 1;
        end
        
        new_order_book = [new_order_book(1:insert_pos-1), order, new_order_book(insert_pos:end)];
    
    else
        if strcmp(order.bs_type, 'buy')
            new_order_book = bid_orders(new_order_book, order);
        else
            new_order_book = ask_orders(new_order_book, order);
        end
    end
end


% used to find whether temporary order price is already in order book
function index = find_price_index(price, order_book)
    prices = arrayfun(@(x) x.price, order_book);
    index = find(prices == price);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% substitution of autocorr function in econometric toolbox：
function [autocorr_val, parcorr_val] = calculate_correlations(data, lags)
    autocorr_val = zeros(1, lags);
    parcorr_val = zeros(1, lags);

    for lag = 1:lags
        shifted_data = data(1:end-lag);
        corr_data = data(lag+1:end);
        
        mean_shifted = mean(shifted_data);
        mean_corr = mean(corr_data);
        
        numerator_auto = sum((shifted_data - mean_shifted) .* (corr_data - mean_corr));
        denominator_auto = sqrt(sum((shifted_data - mean_shifted).^2) * sum((corr_data - mean_corr).^2));
        
        autocorr_val(lag) = numerator_auto / denominator_auto;
        
        for j = 1:lag-1
            shifted_data_j = data(1:end-j);
            corr_data_j = data(j+1:end);
            
            mean_shifted_j = mean(shifted_data_j);
            mean_corr_j = mean(corr_data_j);
            
            numerator_par = sum((shifted_data_j - mean_shifted_j) .* (corr_data_j - mean_corr_j));
            denominator_par = sqrt(sum((shifted_data_j - mean_shifted_j).^2) * sum((corr_data_j - mean_corr_j).^2));
            
            parcorr_val(lag) = parcorr_val(lag) + (numerator_par / denominator_par);
        end
        parcorr_val(lag) = parcorr_val(lag) + autocorr_val(lag);
    end
end



% sort orders by time
function sorted_orders = sort_orders_by_time(orders)
    % check for empty order
    if isempty(orders)
        error('orders is empty');
    end

    % examine the first element in orders: include 'time' or not
    if ~isfield(orders(1), 'time')
        error('orders does not have a ''time'' field');
    end

    % Sort the orders by the 'response_time' field
    [~, order] = sort([orders.response_time]);
    %[~, order] = sort([orders.time]); !!
    sorted_orders = orders(order);
end


% bid order book
% Sort the bid order book by lowest to highest price
% for execution it will search backwards
function bid_order_book = bid_orders(old_book, new_order)

    % !!
    if strcmp(new_order.type, 'market')
        % insert market orders at the last in order to get first executions
        bid_order_book = [old_book, new_order];
        return;
    end

    if isempty(old_book)
        bid_order_book = [old_book, new_order];
        return;
    end

    insert_pos = find([old_book.price] > new_order.price, 1, 'first');

    if isempty(insert_pos)
        % if insert position is empty, insert the new order to the end
        bid_order_book = [old_book, new_order];
    else
        bid_order_book = [old_book(1:insert_pos-1), new_order, old_book(insert_pos:end)];
    end
    % Sort the orders by the 'price' field
    %[~, order] = sort([new_book(:).price]);
    %bid_order_book = new_book(order);
end


% ask order book
% Sort the ask order book by highest to lowest price
% for execution it will search backwards
function ask_order_book = ask_orders(old_book, new_order)
    
    % !!    
    if strcmp(new_order.type, 'market')
        ask_order_book = [old_book, new_order];
        return;
    end

    if isempty(old_book)
        ask_order_book = [old_book, new_order];
        return;
    end

    insert_pos = find([old_book.price] < new_order.price, 1, 'first');    

    if isempty(insert_pos)
        ask_order_book = [old_book, new_order];
    else
        ask_order_book = [old_book(1:insert_pos-1), new_order, old_book(insert_pos:end)];
    end
    % Sort the orders by the 'price' field
    %[~, order] = sort([new_book.price],'descend');
    %ask_order_book = new_book(order);
end


% call auction (execute orders in bid & ask)
% perform call auction to decide opening price and order book
function [new_bid_order_book, new_ask_order_book, exchange_price, investors] = call_auction(bid_order_book, ask_order_book, current_price, investors)
    
    %executed_orders = struct('type', {}, 'price', {}, 'quantity', {}, 'time', {}, 'response_time', {}, 'bs_type', {}, 'agent', {});
    new_ask_order_book = ask_order_book;
    new_bid_order_book = bid_order_book;
    exchange_price = current_price;

    % execute orders
    for i = length(new_bid_order_book):-1:1
        for j = length(new_ask_order_book):-1:1
            if is_matching(new_bid_order_book(i), new_ask_order_book(j))

                exchanged_temp_order = new_bid_order_book(i);
                opposite_order = new_ask_order_book(j);
                exchanged_oppo_order = opposite_order;
                
                % exchange price is best selling price
                exchange_price = opposite_order.price;
                exchanged_temp_order.price = exchange_price;

                % if bid order has smaller or equal quantity
                if abs(new_bid_order_book(i).quantity) <= abs(opposite_order.quantity)
                    exchanged_oppo_order.quantity = - new_bid_order_book(i).quantity;
                    
                    if abs(new_bid_order_book(i).quantity) == abs(opposite_order.quantity)
                        new_ask_order_book(j) = [];
                        % ?
                        % opening_price = (new_bid_order_book(i) + new_ask_order_book(j))/2;
                    else
                        new_ask_order_book(j).quantity = opposite_order.quantity + new_bid_order_book(i).quantity;
                    end
                    new_bid_order_book(i) = [];
                   
                    % @
                    investors(exchanged_oppo_order.agent) = investors(exchanged_oppo_order.agent).record_update(exchanged_oppo_order);
                    investors(exchanged_temp_order.agent) = investors(exchanged_temp_order.agent).record_update(exchanged_temp_order);
                    
                    %executed_orders = [executed_orders, exchanged_temp_order, exchanged_oppo_order];
                    break

                
                % if bid order has larger quantity
                elseif abs(new_bid_order_book(i).quantity) > abs(opposite_order.quantity)
                    exchanged_temp_order.quantity = -opposite_order.quantity;
                    new_bid_order_book(i).quantity = new_bid_order_book(i).quantity - exchanged_temp_order.quantity;
                    new_ask_order_book(j) = [];
                    %executed_orders = [executed_orders, exchanged_temp_order, exchanged_oppo_order];
                    
                    investors(exchanged_oppo_order.agent) = investors(exchanged_oppo_order.agent).record_update(exchanged_oppo_order);
                    investors(exchanged_temp_order.agent) = investors(exchanged_temp_order.agent).record_update(exchanged_temp_order);
                end
            end
        end
    end
end

