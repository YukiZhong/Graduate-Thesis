% taken from JHPS_CPS
function [simulated_rho, xi, f] = estimate_rho(target_population, xi, f)
    
    %target_population = 300;    

    if all(isnan(xi))
        disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    
    %%%%%%
    %data6 = readtable('2010Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'FH:FO');
    %data8 = readtable('2010Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DT:EA');
    %data_wealth = readtable('2010Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'KS:KT');
    % 
    %data6 = readtable('2011Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'CL:CT'); 
    %data8 = readtable('2011Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'BU:CB');
    %data_wealth = readtable('2011Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'LK:LL');

    data6 = readtable('2012Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'BN:BV');
    data8 = readtable('2012Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'BF:BM');
    data_wealth = readtable('2012Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'NE:NF');

    data6 = readtable('2013Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DU:EC');
    data8 = readtable('2013Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DM:DT');
    data_wealth = readtable('2013Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'LI:LJ');

    data6 = readtable('2016Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DR:DZ');
    data8 = readtable('2016Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'EB:EI');
    data_wealth = readtable('2016Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'GY:GZ');

    data6 = readtable('2017Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'CS:DA');
    data8 = readtable('2017Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DC:DJ');
    data_wealth = readtable('2017Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'GH:GI');

    data6 = readtable('2018Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DD:DL');
    data8 = readtable('2018Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DN:DU');
    data_wealth = readtable('2018Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'HB:HC');

    data6 = readtable('2021Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'CW:DE');
    data8 = readtable('2021Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DG:DN');
    data_wealth = readtable('2021Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'HV:HW');
    
    % %2022
    % data6 = readtable('2022Data_JAPAN_new.csv', 'ReadVariableNames', false, 'Range', 'AH:AP');
    % data8 = readtable('2022Data_JAPAN_new.csv', 'ReadVariableNames', false, 'Range', 'AR:AY');
    % data_wealth = readtable('2022Data_JAPAN_new.csv', 'ReadVariableNames', false, 'Range', 'CS:CT');
    % % 
    % % % read answers for A6 and A8
    % data6 = readtable('2023Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'AS:BA');
    % data8 = readtable('2023Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'BC:BJ');
    % data_wealth = readtable('2023Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DL:DM');
    % 
    %data6 = readtable('2024Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'AS:BA');
    %data8 = readtable('2024Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'BC:BJ');
    %data_wealth = readtable('2024Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DN:DO');


    answers8 = table2array(data8);
    answers8(answers8 == 9) = 0;
    answers6 = table2array(data6);
    answers6(answers6 == 9) = 0;
    
    answers_wealth = table2array(data_wealth);
    monthly_wealth = zeros(size(answers_wealth, 1), 1);
    daily_wealth = zeros(size(answers_wealth, 1), 1);

    for i = 1:size(answers_wealth, 1)
        if not(answers_wealth(i, 1) == 999)
            daily_wealth(i) = answers_wealth(i, 1)*10000/20;
            monthly_wealth(i) = answers_wealth(i, 1)*10000;
        elseif not(answers_wealth(i, 2) == 99999)
            daily_wealth(i) = answers_wealth(i, 2)*8;
            monthly_wealth(i) = answers_wealth(i, 2)*8*20;
        else
            daily_wealth(i) = 0;
            monthly_wealth(i) = 0;
        end
    end
    
    %disp(daily_wealth);

    % Find the transition point from Option A to Option B
    X8 = [10,2000,4000,8000,15000,25000,35000,50000];
    X6 = [1000,5000,10000,15000,20000,30000,40000,45000,50000];

    % 2010
    %X8 = [100,500,1000,1500,2000,3000,6000,9000];
    %X6 = [200,500,1000,2000,4000,7000,10000,15000];

    % 2011
    %X8 = [10,2000,4000,8000,15000,25000,35000,50000];
    %X6 = [1000,5000,15000,30000,37500,40000,42500,45000,50000];


    % Define intial values for the MLE
    lb = 0;
    ub = 1;
    initial_params = 0.88;

    options = optimoptions('fmincon', 'TolFun', 1e-6,'TolX', 1e-6, 'Algorithm', 'interior-point', 'MaxIter', 200);
    
    
    % each participant
    
    X8s = [];
    X6s = [];
    for j = 1:size(answers6, 1)

        D6 = answers6(j, :)-1;
        D8 = answers8(j, :)-1;
        W = daily_wealth(j);

        % transition point
        % 0 = buy; 1 = don't buy; -1 = NaN;
        T8 = find(D8==0, 1,"last");
        if ~isempty(T8)
            X8s = [X8s, X8(T8)];
        end

        T6 = find(D6==0, 1,"last");
        if ~isempty(T6)
            X6s = [X6s, X6(T6)];
        end
        


        % call fmincon to optimize

        likelihood_function = @(params) calculate_likelihood(params(1), W, X8, X6, D6, D8);

        [rho_values(j), fval, exitflag, output] = fmincon(likelihood_function, initial_params, ...
            [], [], [], [], lb, ub, [], options);

    end

    disp("=========================");

    counts_X8 = zeros(size(X8));
    for i = 1:length(X8)
        counts_X8(i)= sum(X8s == X8(i));
    end

    counts_X6 = zeros(size(X6));
    for i = 1:length(X6)
        counts_X6(i)= sum(X6s == X6(i));
    end


    % unique_X8 = unique(X8);
    % unique_X6 = unique(X6);
    % 
    % 
    % filtered_counts_X8 = zeros(size(unique_X8));
    % for i = 1:length(unique_X8)
    %     idx = find(counts_X8 == unique_X8(i), 1);
    %     filtered_counts_X8(i) = counts_X8(idx);
    % end
    % 
    % filtered_counts_X6 = zeros(size(unique_X6));
    % for i = 1:length(unique_X6)
    %     idx = find(counts_X6 == unique_X6(i), 1);
    %     filtered_counts_X6(i) = counts_X6(idx);
    % end
    % 
    % 
    % disp(filtered_counts_X8);

    figure;
    subplot(2,1,1);
    bar(X8, counts_X8);
    % histogram(X8s);
    % xlim([0, 50000]);
    % xticks(X8);
    xtickangle(45);
    ylim([0, 400]);
    title('Distribution of Price Transition Point for Lottery Ticket');
    xlabel('Price Level');
    ylabel('Population');
    hold on

    subplot(2,1,2);
    bar(X6, counts_X6);
    % histogram(X6s);
    % xlim([0, 50000]);
    % xticks(X6);
    xtickangle(45);
    ylim([0, 400]);
    title('Distribution of Price Transition Point for Insurance');
    xlabel('Price Level');
    ylabel('Population');
    hold off
    
    % take out outliers
    rho_values(rho_values>1) = 1;
    rho_hat = rho_values(rho_values>0.1);
    mean_r = mean(rho_hat);
    std_r = std(rho_hat);
    disp(mean_r);
    disp(std_r);

    % use ks density
    disp(max(rho_hat));
    [f, xi] = ksdensity(rho_hat);

    % figure;
    % plot(xi, f, 'LineWidth', 2);
    % title('KDE of \rho');
    % xlabel('\rho');
    % ylabel('Density');
    % grid on;

    figure;
    plot(xi, f, 'LineWidth', 2);
    title('KDE of \rho');
    xlabel('\rho');
    xlim([0.5, 1.1]);
    ylabel('Density');
    ylim([0, 20]);
    grid on;

    % generate from KDE
    simulated_rho = randsample(xi, target_population, true, f);
    
    %figure;
    %histogram(simulated_rho, 30, 'Normalization', 'pdf');
    % title('Simulated \rho Distribution');
    % xlabel('\rho');
    % ylabel('Density');
    % grid on;

    disp("~~~~~~~~~~~~~~~~~~~~~~~~~");
    figure;
    histogram(simulated_rho, 30,'Normalization', 'count');
    xlim([0.6, 1]);
    ylim([0, 60]);
    x_limits = xlim;
    y_limits = ylim;
    text(x_limits(2)-0.05*diff(x_limits), y_limits(2)-0.1*diff(y_limits), ...
            {['Mean = ', num2str(mean_r)],['Std = ', num2str(std_r)]}, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12);
    %title('\rho Distribution of agents');
    xlabel('\rho');
    ylabel('frequency');
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',14)

    end

    % generate from KDE
    simulated_rho = randsample(xi, target_population, true, f);

    %disp(mean(simulated_rho));
    %disp(std(simulated_rho));

    %if exist('simulated_parameters_2021.csv', 'file') == 2
        %fileID = fopen(csvFile, 'w');
        %fclose(fileID);
    %else
        %writematrix(simulated_rho', 'simulated_parameters_2021.csv');
    %end
    

    % mle

    function U = crra(x, rho)
        if rho == 1
            U = ln(x);
        else
            U = (x^(1 - rho) - 1)/(1 - rho);
        end
    end

    function w = prospect_weight(p, gamma)
        w = p^gamma/(p^gamma + (1 - p)^gamma)^(1/gamma);
        %w = exp(-(-log(p))^gamma);
    end

    function logL = calculate_likelihood(rho, W, X8, X6, D6, D8)
        logL = 0;

        if W ~= 0
            % lottery
            %w_p = prospect_weight(0.5, gamma);
            
            w_p = 0.5;
            prize = 100000;

            % ~2010
            %prize = 20000;

            for k = 1:length(D8)
                if D8(k) >= 0 && W > X8(k)

                    %disp([X8(k), "X8 calculated"]);
            
                    % define utility funcion
                    EU_lottery = w_p * crra(W + prize - X8(k), rho) + (1-w_p) * crra(W - X8(k), rho);
                    EU_nobuylot = crra(W,rho);
    
                    % define probability function
                    P_lottery = EU_lottery / (EU_lottery + EU_nobuylot);
                    %P_lottery = exp(EU_lottery) / (1 + exp(EU_lottery));
                    %%<-- exceeds calculation capacity

                    L_each_lot = P_lottery^(1-D8(k)) * (1 - P_lottery)^D8(k);
                    
                    
                    %disp(["L_each_lot = ", num2str(L_each_lot(0.5))]);

                    %f_L_lottery = @(rho) L_each_lot(rho) * f_L_lottery(rho);
                    %f_L_lottery = @(rho) L_each_lot(rho) * f_L_lottery(rho);

                    %lottery_funcs{end+1} = L_each_lot;

                    logL = logL + (log(L_each_lot));
                end
            end

            % insurance
            for k = 1:length(D6)

                if D6(k) >= 0 && W > X6(k)

                    EU_insurance = crra(W - X6(k), rho);
                    EU_nobuyins = w_p * crra(max((W - prize), 0),rho) + ...
                        (1-w_p) * crra(W, rho);
    
                    P_insurance = EU_insurance / (EU_insurance + EU_nobuyins);
                    %P_insurance = exp(EU_insurance) / (1 + exp(EU_insurance));
                    %%<-- exceeds calculation capacity
                    
                    L_each_ins = P_insurance^(1-D6(k)) * (1 - P_insurance)^D6(k);
                    
                    %f_L_insurance = @(rho) L_each_ins(rho) * f_L_insurance(rho);
                    %f_L_insurance = @(rho) L_each_ins(rho) * f_L_insurance(rho);

                    %insurance_funcs{end+1} = L_each_ins;
                    logL = logL + (log(L_each_ins));
                end               
            end

        end

        logL = logL;
    end

end