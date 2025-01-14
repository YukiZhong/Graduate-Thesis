function Result = NormFitDemo(Data)
    % error check
    narginchk(1, 1);
    
    Result = struct;
    
    [muhat,sigmahat] = normfit(Data);
    Result.mu = muhat;
    Result.sigma = sigmahat;
    
    % ??kurtosis??
    k = kurtosis(Data);
    Result.kurtosis = k;
    
    % ??skewness??
    s = skewness(Data);
    Result.skewness = s;

    % standard error
    e = std(Data) / sqrt(length(Data));
    Result.stde = e;

    %if nargout == 0
    
    % ------------ creating figure ------------------
    H = figure;
    histfit(Data);
    
    %axis tight;
    % get current axes positions
    ax = gca;
    axPos = ax.Position;
    
    % set string
    str = cell(4,1);
    str{1,1} = ['mu = ',num2str(Result.mu)];
    str{2,1} = ['sigma = ',num2str(Result.sigma)];
    
    % kurtosis
    if Result.kurtosis > 3
        %str{3,1} = ['kurtosis:',num2str(Result.kurtosis),'>3'];
        exc = Result.kurtosis;
        str{3,1} = ['kurtosis = ',num2str(exc), ' (leptokurtic)'];
    elseif Result.kurtosis < 3
        str{3,1} = ['kurtosis = ',num2str(Result.kurtosis),' (platykurtic)'];
    elseif Result.kurtosis == 3
        str{3,1} = 'kurtosis = 3 (mesokurtic)';
    end

    % standard error
    str{4,1} = ['standard error = ', num2str(Result.stde)];
    
    % skewness
    if Result.skewness > 0
        str{5,1} = ['skewness = ',num2str(Result.skewness),' (positive skew)'];
    
    elseif Result.skewness < 0
        str{5,1} = ['skewness = ',num2str(Result.skewness),' (negative skew)'];
    
    elseif Result.skewness == 0
        str{5,1} = 'skewness = 0 (no skew)';
    end
    
    %hold on
    %set(h,'String', str, 'FontSize', 10, 'HorizontalAlignment', 'right', 'BackgroundColor', 'yellow');
    % add text (annotation) on the upper right corner
    annotation('textbox', [axPos(1) + axPos(3) - 0.1, axPos(2) + axPos(4) - 0.1, 0.1, 0.1], ...
               'String', str, 'FitBoxToText', 'on', 'HorizontalAlignment', 'right', 'BackgroundColor', 'yellow', 'FaceAlpha', 0.5);
    %hold off
    %alpha(.5);

end