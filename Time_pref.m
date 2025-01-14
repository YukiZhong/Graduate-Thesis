
function [simulation_results, num_people_each_option, t2_list, sample_size] = Time_pref(target_population, num_people_each_option, t2_list, sample_size)
    %target_population = 300;  % !

    %%%%%%%%%%%%%%%%%% risk-free rate from jgb

    if isnan(num_people_each_option)
    opts = detectImportOptions('jgbcme_all.csv');
    preview('jgbcme_all.csv',opts)

    opts.SelectedVariableNames = {'InterestRate', 'Var11'}; % date and 10Y
    jgb_data = readtable('jgbcme_all.csv', opts);
    jgb_data.Properties.VariableNames = {'Date', 'Value'};
    jgb_data.Date = datetime(jgb_data.Date, 'InputFormat', 'yyyy/MM/dd'); 
    jgb_data.Value = str2double(jgb_data.Value);
    %disp(jgb_data);

    jgb_since_2018 = jgb_data(11178:end, :);  % since 2018
    jgb_since_2010 = jgb_data(9215:12643, :);  % 2010-2024
    jgb_since_2004 = jgb_data(8234:9459, :);  % start of 2006 -2010

    %%%%%%%%%%%%%%%%%%%%
    jgb_subset = jgb_since_2010;

    %opts.VariableTypes = {'datetime', 'double'}; 
    %opts = setvaropts(opts, 'Date', 'InputFormat', 'yyyy-MM-dd');

    length_data = height(jgb_subset)+2;
    disp(length_data);
    
    step_size = 1;
    x = jgb_subset.Date(1:step_size:end);
    y = jgb_subset.Value(1:step_size:end);
    

    figure;
    plot(x, y, "-");   % "-o"
    xlabel('Date');
    ylabel('Value');
    title('Risk-free interest rate from 2010 to 2024');
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',14)


    %%%%%%%%%%%%%%%%%%%%
    %jgb_subset = jgb_since_2018;

    yearly_2007 = round(mean(jgb_data(8482:8726, :).Value)/100, 8);
    yearly_2008 = round(mean(jgb_data(9727:8971, :).Value)/100, 8);
    yearly_2009 = round(mean(jgb_data(8972:9214, :).Value)/100, 8);

    yearly_2011 = round(mean(jgb_data(9460:9704, :).Value)/100, 8);


    yearly_2012 = round(mean(jgb_data(9705:9952, :).Value)/100, 8);
    yearly_2013 = round(mean(jgb_data(9953:10197, :).Value)/100, 8);
    
    
    %yearly_2018 = round(mean(jgb_subset(1:245, :).Value)/100, 8);

    yearly_2016 = round(mean(jgb_data(10686:10930, :).Value)/100, 8);
    yearly_2017 = round(mean(jgb_data(10931:11177, :).Value)/100, 8);

    yearly_2018 = round(mean(jgb_data(11178:11422, :).Value)/100, 8);
    %disp(jgb_subset(1:245, :));
    %disp(yearly_2018);  %    -0.095657%
    
    %yearly_2019 = round(mean(jgb_subset(246:486, :).Value)/100, 8);
    yearly_2019 = round(mean(jgb_data(11423:11663, :).Value)/100, 8);
    %disp(jgb_subset(245:486, :));
    %disp(yearly_2019);  %    -0.095657%

    %yearly_2020 = round(mean(jgb_subset(487:729, :).Value)/100, 8);
    yearly_2020 = round(mean(jgb_data(11664:11906, :).Value)/100, 8);
    %disp(jgb_subset(487:729, :));
    disp(yearly_2020);  %    -0.010417

    yearly_2021 = round(mean(jgb_data(11907:12151, :).Value)/100, 8);
    %disp(jgb_subset(487:729, :));
    disp(yearly_2021);  %    -0.010417

    yearly_2022 = round(mean(jgb_data(12152:12395, :).Value)/100, 8);
    %disp(jgb_subset(487:729, :));
    disp(yearly_2022);  %    -0.010417

    %yearly_2023 = round(mean(jgb_subset(973:1217, :).Value)/100, 8);
    %yearly_2023 = round(mean(jgb_data(12396:12641, :).Value)/100, 8);
    %disp(yearly_2023);  %     0.57092% 

    %%%%%%%%%%%%%%%%% uncertainty (after 13 months)
    % data from 2009~2021
    %data = readtable('JHPS2011data_ver8.0.csv', 'ReadVariableNames', false, 'Range', 'IZ:IZ');
    data = readtable('JHPS2012data_ver8.0.csv', 'ReadVariableNames', false, 'Range', 'JV:JV');
    data = readtable('JHPS2013data_ver8.0.csv', 'ReadVariableNames', false, 'Range', 'JP:JP');
    data = readtable('JHPS2016data_ver7.0.csv', 'ReadVariableNames', false, 'Range', 'NE:NE');
    data = readtable('JHPS2017data_ver6.0.csv', 'ReadVariableNames', false, 'Range', 'ACS:ACS');
    data = readtable('JHPS2018data_ver5.0.csv', 'ReadVariableNames', false, 'Range', 'AAG:AAG');
    % %data = readtable('JHPS2019data_ver4.0.csv', 'ReadVariableNames', false, 'Range', 'ALH:ALH');
    % %data = readtable('JHPS2020data_ver3.0.csv', 'ReadVariableNames', false, 'Range', 'ADQ:ADQ');
    data = readtable('JHPS2021data_ver2.0.csv', 'ReadVariableNames', false, 'Range', 'ADO:ADO');
    %data = readtable('JHPS2022data_ver2.0.csv', 'ReadVariableNames', false, 'Range', 'ADO:ADO');
    % % 
    answers = table2array(data);
    sample_size = length(answers);
    disp(sample_size); 

    figure;
    histogram(answers, 'BinEdges', 0.5:1:8.5, 'Normalization', 'probability');
    xlabel('Options');
    ylabel('Probability');
    title('2021');
    ylim([0, 0.4]);
    set(findall(gcf,'-property','FontSize'),'FontSize',14)

    %10000*(1+r1)^(t1) == 10000*(1+r2)^(t2)

    t1 = 1;
    r1_yearly = yearly_2021;   % ! yearly average



    %r1 = (1+r1_yearly)^(1/12)-1;
    r1 = r1_yearly/12;

    options = [1, 2, 3, 4, 5, 6, 7, 8];
    r2_yearly = [-0.05, 0, 0.02, 0.04, 0.06, 0.10, 0.2, 0.4];
    %r2_monthly = (1+r2_yearly).^(1/12)-1;
    r2_monthly = r2_yearly./12;

    t2_list = [];


    % population for every item
    num_people_each_option = histcounts(answers, [options, inf]);

    % for i = 1:length(options)
    %     if r2_yearly(i) == 0
    %         %random
    %         t2 = NaN;
    %     elseif sign(r2_yearly(i)) ~= sign(r1)
    %         % assume r1 = -10%
    %         t2 = t1*log(1-r1*3)/log(1+r2_monthly(i));
    %         %t2 = t1*log(1-r1*(rand+0.2)*10)/log(1+r2_monthly(i));
    %     else
    %         t2 = t1*log(1+r1)/log(1+r2_monthly(i));
    %     end
    % 
    %     % convert to days
    %     t2 = t2 * 30;
    % 
    %     t2_list = [t2_list, t2];
    % end

    for i = 1:length(options)
        if r2_yearly(i) == 0
            %random
            t2 = NaN;
        elseif r2_yearly(i) < 0 && r1 < 0
            t2 = t1*r1/r2_monthly(i);
        elseif r2_yearly(i) > 0 && r1 < 0
            t2 = -t1*r1/log(1+r2_monthly(i));

        elseif r2_yearly(i) < 0 && r1 > 0
            t2 = t1*log(1+r1)/-r2_monthly(i);
        else
            t2 = t1*log(1+r1)/log(1+r2_monthly(i));
        end

        % convert to days
        t2 = t2 * 30;

        t2_list = [t2_list, t2];
    end

    disp(t2_list);
    
    end
    
    simulation_results = [];
    simulated_sample = round(num_people_each_option * target_population/sample_size);
    
    if sum(simulated_sample) ~= target_population
        difference = target_population - sum(simulated_sample);
        simulated_sample(end) = simulated_sample(end) + difference;
    end

    
    % distribute t2 to all population
    scale = 1 / t2_list(end);
    scale = 20;
    disp(scale);

    % random generation
    upper_bound = t2_list(3) * scale;
    lower_bound = 1;
    mu = mean(t2_list(3:end)) * scale;
    sigma = std(t2_list(3:end)) * scale;

    for i = 1:length(t2_list)
        if i == 2
            t2_replicated = round(max(normrnd(mu, sigma, [simulated_sample(i), 1]), 1));
        else
            t2_replicated = repmat(round(t2_list(i)*scale), simulated_sample(i), 1);
        end
        simulation_results = [simulation_results; t2_replicated];
    end

    mean_r = mean(simulation_results);
    std_r = std(simulation_results);
    

    figure;
    histogram(simulation_results,10,'Normalization', 'count');
    xlim([0, 100]);
    ylim([0, 100]);
    x_limits = xlim;
    y_limits = ylim;
    text(x_limits(2)-0.05*diff(x_limits), y_limits(2)-0.1*diff(y_limits), ...
            {['Mean = ', num2str(mean_r)],['Std = ', num2str(std_r)]}, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 14);
    xlabel('Time horizon');
    ylabel('Population');
    %title('Simulated time horizon distribution among agents\n');
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',14)

    %if exist('simulated_parameters_2021.csv', 'file') == 2
        %writematrix(simulation_results, 'simulated_parameters_2021.csv', 'WriteMode', 'append');
    %else
        %writematrix(simulation_results, 'simulated_parameters_2021.csv');
    %end
   
end