% taken from JHPS_CPS
function rho_all()
    
    data6_2012 = readtable('2012Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'BN:BV');
    data8_2012 = readtable('2012Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'BF:BM');
    data_wealth_2012 = readtable('2012Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'NE:NF');
    % % 
    data6_2013 = readtable('2013Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DU:EC');
    data8_2013 = readtable('2013Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DM:DT');
    data_wealth_2013 = readtable('2013Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'LI:LJ');
    % % 
    data6_2016 = readtable('2016Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DR:DZ');
    data8_2016 = readtable('2016Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'EB:EI');
    data_wealth_2016 = readtable('2016Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'GY:GZ');
    % % 
    data6_2017 = readtable('2017Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'CS:DA');
    data8_2017 = readtable('2017Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DC:DJ');
    data_wealth_2017 = readtable('2017Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'GH:GI');
    % 
    % %2018
    data6_2018 = readtable('2018Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DD:DL');
    data8_2018 = readtable('2018Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DN:DU');
    data_wealth_2018 = readtable('2018Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'HB:HC');
    % 
    % %2021
    data6_2021 = readtable('2021Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'CW:DE');
    data8_2021 = readtable('2021Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DG:DN');
    data_wealth_2021 = readtable('2021Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'HV:HW');
    % 
    % %2022
    data6_2022 = readtable('2022Data_JAPAN_new.csv', 'ReadVariableNames', false, 'Range', 'AH:AP');
    data8_2022 = readtable('2022Data_JAPAN_new.csv', 'ReadVariableNames', false, 'Range', 'AR:AY');
    data_wealth_2022 = readtable('2022Data_JAPAN_new.csv', 'ReadVariableNames', false, 'Range', 'CS:CT');
    % % 
    % % % read answers for A6 and A8
    data6_2023 = readtable('2023Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'AS:BA');
    data8_2023 = readtable('2023Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'BC:BJ');
    data_wealth_2023 = readtable('2023Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DL:DM');
    % 
    data6_2024 = readtable('2024Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'AS:BA');
    data8_2024 = readtable('2024Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'BC:BJ');
    data_wealth_2024 = readtable('2024Data_JAPAN.csv', 'ReadVariableNames', false, 'Range', 'DN:DO');

    
    allDataTable6 = table();
    allDataTable8 = table();
    allDataTable_wealth = table();

    years = {'2012', '2013', '2016', '2017', '2018', '2021', '2022'};

    for j = 1:length(years)
        data6 = eval(['data6_', years{j}]);
        data8 = eval(['data8_', years{j}]); 
        data_wealth = eval(['data_wealth_', years{j}]);

        answers8 = table2array(data8);
        answers8(answers8 == 9) = 0;
        answers6 = table2array(data6);
        answers6(answers6 == 9) = 0;

        answers_wealth = table2array(data_wealth);

        daily_wealth = [];
        for i = 1:size(answers_wealth, 1)
            if answers_wealth(i, 1) ~= 999 && answers_wealth(i, 1) ~= 0
                daily_wealth = [daily_wealth, answers_wealth(i, 1)*10000/20];
            elseif answers_wealth(i, 2) ~= 99999 && answers_wealth(i, 2) ~= 0
                daily_wealth = [daily_wealth, answers_wealth(i, 2)*8];
            end
        end


        % find transition point 8 = lottery
        X8 = [10,2000,4000,8000,15000,25000,35000,50000];
        X6 = [1000,5000,10000,15000,20000,30000,40000,45000,50000];
    
        X8s = [];
        X6s = [];
        for k = 1:size(answers6, 1)
    
            D6 = answers6(k, :)-1;
            D8 = answers8(k, :)-1;
    
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
        end


        
        currentyear = str2double(years{j});
        disp([repmat(currentyear, size(X6s))', X6s']);
        yearTable6 = array2table([repmat(currentyear, size(X6s))', X6s'], 'VariableNames', {'Year', 'ConversionPoint'}); % 创建表格
        allDataTable6 = [allDataTable6; yearTable6]; 

        yearTable8 = array2table([repmat(currentyear, size(X8s))', X8s'], 'VariableNames', {'Year', 'ConversionPoint'}); % 创建表格
        allDataTable8 = [allDataTable8; yearTable8]; 

        yearTable_wealth = array2table([repmat(currentyear, size(daily_wealth))', daily_wealth'], 'VariableNames', {'Year', 'Wealth'}); % 创建表格
        allDataTable_wealth = [allDataTable_wealth; yearTable_wealth];
    end
    
    %disp(allDataTable);
    

    % 2010
    %X8 = [100,500,1000,1500,2000,3000,6000,9000];
    %X6 = [200,500,1000,2000,4000,7000,10000,15000];

    % 2011
    %X8 = [10,2000,4000,8000,15000,25000,35000,50000];
    %X6 = [1000,5000,15000,30000,37500,40000,42500,45000,50000];

    
    disp("=========================");
    % boxplots

    data = allDataTable6{:, 'ConversionPoint'};
    groups = allDataTable6{:, 'Year'};
    
    figure;
    boxplot(data, groups, 'Labels', unique(groups));

    xlabel('Year');
    ylabel('Conversion Point of Purchasing Insurance');
    set(findall(gcf,'-property','FontSize'),'FontSize',14)

    % =========================

    data = allDataTable8{:, 'ConversionPoint'};
    groups = allDataTable8{:, 'Year'};
    
    figure;
    boxplot(data, groups, 'Labels', unique(groups));

    xlabel('Year');
    ylabel('Conversion Point of Purchasing Lottery Ticket');
    set(findall(gcf,'-property','FontSize'),'FontSize',14)

    % =========================

    data = allDataTable_wealth{:, 'Wealth'};
    groups = allDataTable_wealth{:, 'Year'};
    
    figure;
    boxplot(data, groups, 'Labels', unique(groups));

    xlabel('Year');
    ylabel('Daily Income of Participants');
    ylim([0, 55000]);
    set(findall(gcf,'-property','FontSize'),'FontSize',14)


end