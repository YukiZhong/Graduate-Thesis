classdef me_agents
    properties
        id;
        response_time_mu;
        investor_type;
        initial_wealth;
        time_steps;

        n_F;
        n_P;
        n_N;
        tau;    %inv_memory

        % values that change for each t
        holding = 0;    %inv_stock
        cash_balance;   %inv_money
        response_time;

        p_bar;
        predicted_price;
        position;
        order_type;
        quote_price_t;
        quantity_t;
        bs_type;

        % fixed parameters
        %r = 0.002; % risk-free interest rate
        r = 0.0012;

        %delta = 0.01; % CARA,!change to 0.01
        delta = 0;
        tick_size = 0.01;

        %tau_0 = 10; % Average time level (common for all agents)
        tau_0;
        alpha_0;

        max_quantity = 300;   % !!!
        absolute_risk_aversion_CARA = 1;
        % !!
        gamma_= 0.1;   % Volatility to estimate position ! change to 0.002
        sigma_D = 0.002; % Variance of reaction times
        
        %utility_ts = []; % initialise utility
        %positions = [];
        %profits = [];
        wealth = [];
        %executed = [];
    end
    
    methods
        % object construction
        function obj = me_agents(id, mu, initial_wealth, wealth, type, num_steps, initial_stock_price, tau_0, alpha_0)
            obj.id = id;
            obj.response_time_mu = mu;
            obj.initial_wealth = initial_wealth;
            obj.wealth = wealth;
            obj.investor_type = type;
            obj.time_steps = num_steps;
            %obj.utility_ts = zeros(1, num_steps); % initialise utility
            %obj.positions = zeros(1, num_steps);
            obj.quote_price_t = initial_stock_price;
            obj.cash_balance = initial_wealth;
            obj.tau_0 = tau_0;
            obj.alpha_0 = alpha_0;
        end
        
        
        % reaction time at each t
        function obj = cal_time(obj, t)
            obj.response_time = t + rand;
        end

        
        % tau for memory length
        function obj = calculate_tau(obj, n_F, n_P, n_N)
            obj.n_F = n_F;
            obj.n_P = n_P;
            obj.n_N = n_N;
            obj.tau = max(round(obj.tau_0 * (1 + n_F) / (1 + n_P)), 1);
        end

        
        % average price over history prices
        function obj = calculate_p_bar(obj, stock_prices, t)
            if t == 1
                obj.p_bar = stock_prices(1);
            elseif t <= obj.tau
                obj.p_bar = mean(stock_prices(1:t));
            else
                obj.p_bar = mean(stock_prices(t+1-obj.tau:t));
            end
        end

        
        % predict
        function obj = generate_predicted_price(obj, p_F, epsilon, current_price)
            obj.predicted_price = (obj.n_F * p_F + obj.n_P * obj.p_bar + obj.n_N * exp(epsilon) * current_price) / (obj.n_F + obj.n_P + obj.n_N);
        
            % ensure predicted price is greater than 0
            if obj.predicted_price <= 0
                obj.predicted_price = obj.tick_size;
            end
        end


        % decide quote price first before calculating position
        function obj = quote_price(obj, bid_order_book, ask_order_book, current_price)
            quote_price = current_price;

            % set price for limit order
            if strcmp(obj.order_type, 'limit')
                %if quantity > 0
                    % buy
                    best_sell_price = me_agents.get_best_order_price(ask_order_book, current_price);
                    if best_sell_price < 0
                        error('Best sell price is negative');
                    end
                    %quote_price = max(obj.tick_size, min(round(obj.predicted_price, 2), round(best_sell_price, 2)));
                    
                %elseif quantity < 0
                    % sell
                    best_buy_price = me_agents.get_best_order_price(bid_order_book, current_price);
                    if best_buy_price < 0
                        error('Best buy price is negative');
                    end
                    %quote_price = max(obj.tick_size, max(round(obj.predicted_price, 2), round(best_buy_price, 2)));
                %end
                
                %u = rand(1)*5*obj.tick_size; 
                u = 20*obj.tick_size;
                lower_bound = best_buy_price - u;
                upper_bound = best_sell_price + u;

                % obj.exch_price=round((rand(1)*(n-m)+m)*100)/100;
                quote_price = round((lower_bound + (upper_bound - lower_bound) * rand(1))*100)/100;
            end
            
            obj.quote_price_t = quote_price;
        end



        % calculate positions of investors
        function obj = calculate_position(obj, t, stock_prices)
            
            % Calculate position for CARA investor
            if strcmp(obj.investor_type, 'CARA')
                % Initial average of risk aversion

                %alpha_0 = 0.8;
                alpha = obj.alpha_0 * (1 + obj.n_F) / (1 + obj.n_P);
               
                %if t > 1
                    %previous_price = stock_prices(t-1);
                %else
                    previous_price = stock_prices(t); % or another default value
                %end
                
                % calculate variance of p
                p_var = me_agents.calculate_p_variance(stock_prices, obj.tau, t, obj.gamma_);
                pos = (obj.predicted_price - ((1 + obj.r)^obj.tau) * (previous_price + obj.delta)) / (alpha * p_var);
                %position = (predicted_price - (1 + r)^tau * (stock_prices(t) + gamma_)) / (alpha * p_var);
        
                % if t == 2
                %     disp(stock_prices(1:4));
                %     disp(obj.tau);
                %     disp(p_var);
                % end

            % Calculate position for IARA investor
            elseif strcmp(obj.investor_type, 'IARA')
                % Initial average of risk aversion
                % a is just a parameter
                alpha_1 = 0.02;
                alpha_1 = 0.012;
                alpha_1 = 0.0175;

                a = obj.alpha_0 * (1 + obj.n_F) / (1 + obj.n_P);

                if t > 1
                    %previous_price = stock_prices(t-1);
                    previous_price = stock_prices(t-1);
                    f_t = obj.wealth(t-1); % wealth holding by investor i at time t
                %elseif t > 2
                    %previous_price = stock_prices(t-2);
                    %f_t = obj.wealth(t-1); % wealth holding by investor i at time t
                else
                    previous_price = stock_prices(t); % or another default value
                    f_t = obj.initial_wealth; % wealth holding by investor i at time 1
                end
        
                %lambda = (2*a)/((-2)*a*f_t+1);
                
                % calculate variance of p
                p_var = me_agents.calculate_p_variance(stock_prices, obj.tau, t, obj.gamma_);
        
                %pos = ((1 - 2 * a * f_t) * (obj.predicted_price-obj.quote_price_t) / (2 * a * p_var + 2 * a * (obj.predicted_price)^2 - 2 * a * obj.predicted_price * previous_price));
                pos = ((1 - 2 * a * f_t) * (obj.predicted_price-previous_price)) / (2 * a * p_var + 2 * a * (obj.predicted_price)^2 - 2 * a * obj.predicted_price * previous_price);
                
                %disp('position: ');
                %disp(pos);

                if isnan(pos)
                    disp(obj.id);
                    disp('numerator');
                    disp((1 - 2 * a * f_t) * (obj.predicted_price-obj.quote_price_t));
                    disp('quote_price_t');
                    disp(obj.quote_price_t);
                    disp('denominator');
                    disp(2 * a * p_var + 2 * a * obj.predicted_price^2 - 2 * a * obj.predicted_price * previous_price);
                    disp(p_var);
                    error('pos is nan');
                end
        
            else
                % Calculate position for Random investor
                pos = round(randn * 10, 2);
                % randn has mean 0 and variance 1
            end

            obj.position = pos;

            obj.quantity_t = me_agents.decide_limit_order_params(obj.position, obj.max_quantity, obj.holding);
        end
        
                
        % randomly decide the order type base on market_order_ratio
        function obj = decide_order_type(obj, current_market_order_ratio, call)
            if call
                obj.order_type = 'limit';
            elseif ~call
                if rand < current_market_order_ratio
                    obj.order_type = 'market';
                else
                    obj.order_type = 'limit';
                end
            end
            %fprintf("id: %d; order type: %s\n", obj.id, obj.order_type);
        end
        
        % create new orders
        function new_order = create_orders(obj, t)
            % Create a new order
            if obj.quantity_t > 0
                obj.bs_type = 'buy';
            else
                obj.bs_type = 'sell';
            end
            new_order = struct('type', obj.order_type, 'price', obj.quote_price_t, 'quantity', obj.quantity_t, 'time', t, 'response_time', obj.response_time, 'bs_type', obj.bs_type, 'agent', obj.id);
        end
    

        
        % calculate remaining cash_balance
        function obj = record_update(obj, current_order)
    
            % check for errors
            if current_order.quantity == 0 || isnan(current_order.quantity)
                error('executed_orders(i).quantity is zero or NaN');
            end
    
            if current_order.price == 0
                error('executed_orders(i).price is zero');
            elseif isnan(current_order.price)
                error('executed_orders(i).price is nan');
            end
    
            % spend for buy orders; and income for sell orders
            obj.holding = obj.holding + current_order.quantity;
            obj.cash_balance = obj.cash_balance - current_order.quantity * current_order.price;
        end


        
        % update wealth and utility at the end of each time step
        function obj = calculate_wealth_t(obj, t, current_stock_price)
            
            % Calculate the total wealth for now
            %disp("----------------------");
            %disp(obj.cash_balance);
            %disp(obj.holding);
            %disp(current_stock_price);
            obj.wealth(t) = obj.cash_balance + obj.holding * current_stock_price;
            %obj.profits(t) = obj.wealth(t) - obj.initial_wealth;

            % Calculate the utility of the total wealth
            
            %utility = me_agents.calculate_utility(obj.wealth(t), obj.investor_type, obj.absolute_risk_aversion_CARA, obj.initial_wealth, obj.alpha_0);
            %obj.utility_ts(t) = utility;
        end


        % calculate final profit after all time steps
        function final_wealth = cal_final_profit(obj)
            %disp(obj.wealth);
            final_wealth = obj.wealth(end) - obj.initial_wealth;
        end

    end

    
    
    methods (Static)

        function quantity_t = decide_limit_order_params(position, max_quantity, holding)

            quantity_t = round(position) - holding;
            quantity_t = round(quantity_t, 2);

            if abs(quantity_t) > max_quantity
                quantity_t = max_quantity * sign(quantity_t);
            end
        end


        function best_order_price = get_best_order_price(order_book, current_price)
            % if the opposite side order_book is empty，return current price
            if isempty(order_book)
                best_order_price = current_price;
                return;
            end            

            % determine best order price by order type
            if strcmp({order_book.bs_type}, 'buy')
                best_order_price = max([order_book.price]);
            else
                best_order_price = min([order_book.price]);
            end
        end

        % calculate variance of p for all type of agents
        function p_var = calculate_p_variance(stock_prices, tau, t, gamma_)
            if t < 2 || tau < 2  % ! if t < 2 
                p_var = gamma_; 
                p_var2 = gamma_;
            elseif t <= tau
                p_var = var(stock_prices(1:t));

                log_returns = [];
                for a = 2:t
                    next = log(stock_prices(a)/stock_prices(a-1));
                    log_returns = [log_returns,next];
                end
                p_var2 = var(log_returns);

            else
                p_var = var(stock_prices(t-tau+1:t));
                
                log_returns = [];
                for a = t+2-tau:t
                    next = log(stock_prices(a)/stock_prices(a-1));
                    log_returns = [log_returns,next];
                end
                p_var2 = var(log_returns);

            end
        
            % set a safe guardance
            if p_var < 0.01
                p_var = 0.01;
            end

            if p_var2 < 0.01
                p_var2 = 0.01;
            end
        end


        % calculat utility
        function utility = calculate_utility(wealth, investor_type, absolute_risk_aversion_CARA, initial_wealth_IARA, alpha_0)
            if wealth <= 0
                utility = -inf;
            elseif strcmp(investor_type, 'CARA')
                % Calculate utility for CARA investor
                utility = -exp(-absolute_risk_aversion_CARA * wealth);
            elseif strcmp(investor_type, 'IARA')
                % Calculate utility for IARA investor
                utility = (wealth - initial_wealth_IARA) / (1 - alpha_0);
            else
                % Calculate utility for Random investor
                utility = log(wealth);
            end
        end

    end
end