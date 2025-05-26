clear
clc

data = readtable(['Openloop/Openloop_Step_Tg_P.csv']);
data2 = readtable(['Openloop/Openloop_Step_Tg_N.csv']);
%data3 = readtable(['TempFlowFeedback/Step_Tg_N.csv']);

%% Plot multiple
time_indices = [1000]; 
savepath = false;

plot_index = false;
if plot_index
    plot_all_of_var(data, 'T', 'K', time_indices, savepath);
    plot_all_of_var(data, 'Tg','K', time_indices, savepath);
    plot_all_of_var(data, 'M', 'kg',time_indices, savepath);
    plot_all_of_var(data, 'H', 'kJ',time_indices, savepath);
    plot_all_of_var(data, 'beta','-', time_indices, savepath);
    plot_all_of_var(data, 'Q', 'kW',time_indices, savepath);
    plot_all_of_var(data, 'm', 'kg/s',time_indices, savepath);
    plot_all_of_var(data, 'p', 'bar',time_indices, savepath);
    plot_all_of_var(data, 'mg', 'kg/s', time_indices, savepath);
end


%% Plot Single
var = { 
    {'m_0', 'kg/s'}, {'M_0', 'kg'}, {'T_0', 'K'}, {'mg_0', 'kg/s'}, {'Tg_0', 'K'}, ...
    {'Q_0', 'kW'}, {'beta_0', '-'}, {'p_0', 'bar'}, {'H_0', 'kJ'}, ...
    {'m_9', 'kg/s'}, {'M_9', 'kg'}, {'T_9', 'K'}, {'mg_9', 'kg/s'}, {'Tg_9', 'K'}, ...
    {'Q_9', 'kW'}, {'beta_9', '-'}, {'p_9', 'bar'}, {'H_9', 'kJ'}, ...
    {'m_19', 'kg/s'}, {'M_19', 'kg'}, {'T_19', 'K'}, {'mg_19', 'kg/s'}, {'Tg_19', 'K'}, ...
    {'Q_19', 'kW'}, {'beta_19', '-'}, {'p_19', 'bar'}, {'H_19', 'kJ'}, ...
    {'m_29', 'kg/s'}, {'M_29', 'kg'}, {'T_29', 'K'}, {'mg_29', 'kg/s'}, {'Tg_29', 'K'}, ...
    {'Q_29', 'kW'}, {'beta_29', '-'}, {'p_29', 'bar'}, {'H_29', 'kJ'}, ...
    {'m_36', 'kg/s'}, {'M_36', 'kg'}, {'T_36', 'K'}, {'mg_36', 'kg/s'}, {'Tg_36', 'K'}, ...
    {'Q_36', 'kW'}, {'beta_36', '-'}, {'p_36', 'bar'}, {'H_36', 'kJ'}, ...
};



var2 = {{'m_in', 'kg/s'}, {'m_s', 'kg/s'}, {'Ts', 'K'}, {'TBs', 'K'}, {'pS', 'bar'}, ... 
     {'hS', 'kJ/kg'}, {'beta_s', '-'}, {'rho_s', 'kg/m^3'}, ...
     {'m_T', 'kg/s'}, {'TT', 'K'}, {'TBT', 'K'}, {'pT', 'bar'}, ...
     {'hT', 'kJ/kg'}, {'beta_T', '-'}, {'rho_T', 'kg/m^3'}, {'P', 'kW'}, ...
     {'Ms', 'kg'}, {'Hs', 'kJ'}, {'MT', 'kg'}, {'HT', 'kJ'}, {'dummy_zV', '-'} ...
     {'dummy_zm', '-'}, {'dummy_mg', 'kg/s'}};

var3 = {{'m_U', 'kg/s'}, {'TU', 'K'}, {'TBU', 'K'}, {'hU', 'kJ/kg'}, ...
     {'beta_U', '-'}, {'Tc', 'K'}, {'m_c', 'kg/s'}, ...
     {'hc', 'kJ/kg'}, {'Qc', 'kW'}
    };

var4 = {{'P', 'kW'},{'pS', 'bar'}, {'dummy_zm', '-'}, {'Ts', 'K'}, ...
    {'dummy_zV', '-'}, {'dummy_mg', 'kg/s'}, {'dummy_pS_SP', 'bar'}...
    {'dummy_pS_SP_VPC','bar'}};

var5 = {{'Ts', 'K'}};

%{'dummy_v', '-'},
end_time = 2000;
%plot_propagation(data, 'm', 'kg/s', [0,9,19,29,36], false);
%temperature_gradient(data, savepath, {'T', 'Tg'}, time_indices);
%plot_single_hex(data, var, savepath, end_time);
%plot_single(data, var5, savepath, end_time);
%plot_single(data, var3, savepath, end_time);
%plot_single(data, var4, savepath, end_time);
plot_double(data, data2, var5, savepath, end_time);
%plot_triple(data, data2, data3, var5, savepath, end_time);
%plot_triple(data, data2, data3, var3, savepath, end_time);
%plot_triple(data, data2, data3, var4, savepath, end_time);
%step = [(10.95)*ones(1,200), (10.55)*ones(1,1800)];
%step2 = [(443.338+273.15)*ones(1,200), (443.338+273.15-10)*ones(1,1800)];
%plot_step(step);
%plot_step_two(step, step2);
%plot_two(data, {'beta_16', 'beta_42'}, savepath, end_time)

%plot_cross_beta(data, {'beta_27'}, savepath, end_time)
%N = 37; % Number of beta variables (0 to 36)

% Generate crossVars
%crossVars = arrayfun(@(i) ['beta_', num2str(i)], 0:N-1, 'UniformOutput', false);

% Generate crossingTargets
%crossingTargets = zeros(1,N); %[zeros(1, floor(N/2)), ones(1, ceil(N/2))];

%plot_cross(data, {'M_10', 'kg'}, savepath, end_time, crossVars, crossingTargets)

%rmse_deviation(data, {'Ts', 'K'}, end_time, 682)

%% Plot functions
function f = plot_all_of_var(data, variable_prefix, unit, time_indices, savepath)
    variable_names = data.Properties.VariableNames;
    colorss = {  '#333CF5','#FFB80E','#BC19BF', '#1CB6CC', '#10DC28'};
    % Initialize a figure
    figure;
    hold on;
    
    % Loop through each time index
    for t_idx = 1:length(time_indices)
        t = time_indices(t_idx);
        % Initialize arrays to store indices and values

        indices = [];
        values = [];

        % Loop through each variable and extract the value at the specified time index
        for i = 1:length(variable_names)
            var_name = variable_names{i};
            if startsWith(var_name, [variable_prefix '_'])
                % Extract the index from the variable name
                index = str2double(regexp(var_name, '\d+', 'match'));
                value = data{t, var_name};
                indices = [indices, index];
                values = [values, value];
            end
        end
        
        % Sort the indices and corresponding values
        [indices, sortIdx] = sort(indices);
        values = values(sortIdx);

        color = colorss{mod(t_idx-1, length(colorss)) + 1};
        
        % Plot the data
        plot(indices+1, values, '-o', 'DisplayName', ['t = ', num2str(t),'s'], 'Color',color, 'Linewidth', 0.8);
    end

    % Customize the plot
    xlabel('Index');
    if strcmp(variable_prefix, 'mg')  
           ylabel(['m_g', ' [', unit, ']']);
    elseif strcmp(variable_prefix, 'm')  
           ylim([10,11]);
           ylabel([variable_prefix, ' [', unit,']']);
    elseif strcmp(variable_prefix, 'beta')  
           ylabel(['\beta', ' [', unit,']']);
    elseif strcmp(variable_prefix, 'Tg')  
           ylabel(['T_g', ' [', unit,']']);
    else
        ylabel([variable_prefix, ' [', unit,']']);
    end
    %legend show;
    %legend('Location', 'best');
    xlim([1,37]);

    grid on;
    hold off;
    if savepath
        if variable_prefix == 'm'
            saveas(gcf, [savepath,'massflow','.png']);
        else
            saveas(gcf, [savepath,variable_prefix,'.png']);
        end
    end
    f = true;
end


function f = plot_single_hex(data, var, savepath, end_time)
    numElements = length(var);
    
    i = 1;

    while true
        if i > numElements
            break;
        end

        currentVar = var{i}{1}; % Variable name
        unit = var{i}{2};       % Unit string

        yData = data.(currentVar);
       
        yDataFiltered = yData(1:end_time); 
        timeFiltered = linspace(0, length(yDataFiltered), length(yDataFiltered)); 
        
        % correct for 0 indexing in casadi
        parts = split(currentVar, '_');
        segmentVar = parts{1};
        segmentIndex = str2double(parts{2}) + 1;
        segmentIndex = num2str(segmentIndex);

        fig = figure('Visible', 'on');

        hold on
        plot(timeFiltered, yDataFiltered, 'color', '#333CF5', 'LineWidth', 0.8);
        ylabel([segmentVar, '_{', segmentIndex, '} [', unit, ']']);
        xlabel('Time [s]');
        xlim([0,end_time]);

        %if segmentVar == 'm'  
        %    ylim([10,11]);
        %end

        grid on;
        hold off;

        if savepath
            if segmentVar == 'm'
                saveas(fig, [savepath,'massflow', '_', segmentIndex,'.png']);
            else
                saveas(fig, [savepath, segmentVar, '_', segmentIndex,'.png']);
            end
        end

        %close(fig);
    
        i = i + 1;
    end


    f = true;
end

function f = plot_single(data, var, savepath, end_time)
    numElements = length(var);
    
    i = 1;

    while true
        if i > numElements
            break;
        end

        currentVar = var{i}{1}; % Variable name
        unit = var{i}{2};       % Unit string

        yData = data.(currentVar);
       
        yDataFiltered = yData(1:end_time); 
        timeFiltered = linspace(0, length(yDataFiltered), length(yDataFiltered)); 
        
        % correct for 0 indexing in casadi

        fig = figure('Visible', 'on');

        hold on
        plot(timeFiltered, yDataFiltered, 'color', '#333CF5', 'LineWidth', 0.8);
        %ylabel([currentVar, '[', unit, ']']);
        xlabel('Time [s]');
        xlim([0,end_time]);

        if strcmp(currentVar, 'TU_real')  
           ylabel(['TU_{real}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'm_cold_out')
            ylabel(['m_{cold}^{out}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'T_cold_out')
            ylabel(['T_{cold}^{out}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'm_in')
            ylabel(['m_{p}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hS')
            ylabel(['h_S', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hT')
            ylabel(['h_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'pS')
                ylabel(['p_s', ' [', unit, ']']);
                %yData2 = data.('dummy_pS_SP');
                %yData2Filtered = yData2(1:end_time); 
                %plot(timeFiltered, yData2Filtered, 'LineStyle','--', 'DisplayName','Setpoint', 'LineWidth', 0.8, 'color', '#FFB80E');
                %plot(timeFiltered, 22*ones(1, end_time), 'LineStyle','--', 'Color','black', 'DisplayName','Boundary', 'LineWidth', 0.8);
                %plot(timeFiltered, 24*ones(1, end_time), 'LineStyle','--', 'Color','black', 'LineWidth', 0.8);
                %ylim([21.9,24.1]);
                %legend('p_s', 'Setpoint', 'Boundaries');
                %legend('location', 'best');
        elseif strcmp(currentVar, 'pT')
                ylabel(['p_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'Ts')
                ylabel(['T_s', ' [', unit, ']']);
                %plot(timeFiltered, 682*ones(1, end_time), 'LineStyle','--', 'Color','black', 'DisplayName','Setpoint');
        elseif strcmp(currentVar, 'P')
            ylabel(['P', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TT')
                ylabel(['T_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_zm')
                ylabel(['z_m', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_mg')
                ylabel(['m_g^0', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_zV')
                ylabel(['z_v', ' [', unit, ']']);
        elseif strcmp(currentVar, 'HT')
                ylabel(['H_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hU')
                ylabel(['h_U', ' [', unit, ']']);
        elseif strcmp(currentVar, 'MT')
                ylabel(['M_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'Ms')
                ylabel(['M_s', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TBs')
                ylabel(['T_s^{sat}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TBT')
                ylabel(['T_T^{sat}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TBU')
                ylabel(['T_U^{sat}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TU')
                ylabel(['T_U', ' [', unit, ']']);
        elseif strcmp(currentVar, 'rho_s')
            ylabel(['\rho_s', ' [', unit, ']']);
        elseif strcmp(currentVar, 'rho_T')
            ylabel(['\rho_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_m_in_SP')
            ylabel(['m_{p}^{SP}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hc')
            ylabel(['h_c', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_pS_SP')
            ylabel(['p_s^{SP}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_pS_SP_VPC')
            ylabel(['p_{s,0}^{SP}', ' [', unit, ']']);         
        else 
           ylabel([currentVar, ' [', unit, ']']);
        end

        grid on;
        hold off;

        if savepath
            if currentVar == 'm'
                saveas(fig, [savepath,'massflow', '.png']);
            elseif currentVar == 'P'
                saveas(fig, [savepath,'Power', '.png']);
            elseif strcmp(currentVar, 'hS')
                saveas(fig, [savepath,'h_S', '.png']);
            elseif strcmp(currentVar, 'hT')
                saveas(fig, [savepath,'h_T', '.png']);
            elseif strcmp(currentVar, 'hc')
                saveas(fig, [savepath,'h_c', '.png']);
            else
                saveas(fig, [savepath, currentVar, '.png']);
            end
        end

        %close(fig);
    
        i = i + 1;
    end


    f = true;
end

function f = plot_step(data)
    
    yData = data;
    time = linspace(0, 2000, 2000); 
    
    figure;
    hold on
    plot(time, yData, 'Color', 'black', 'LineWidth', 1);
    ylabel('m_p^{SP} [kg/s]');
    xlabel('Time [s]');
    xlim([0, 2000]);
    ylim([10.5,11]);
    grid on;
    hold off;

    saveas(gcf, 'TempFlowFeedback/ConstP/step.png');

    f = true;
end

function f = plot_step_two(data, data2)
    
    yData = data;
    yData2 = data2;
    time = linspace(0, 2000, 2000); 
    
    figure;
    hold on
    plot(time, yData, '-', 'color', '#333CF5', 'LineWidth',0.8);
    plot(time, yData2, '--', 'color', '#FFB80E', 'LineWidth',0.8);
    ylabel('T_g^0 [K]');
    xlabel('Time [s]');
    xlim([0, 2000]);
    %ylim([22.9, 23.1]);
    legend('T_g^0 + 10K', 'T_g^0 - 10K');
    legend('Location', 'best');
    grid on;
    hold off;

    saveas(gcf, 'Openloop/Nonlinearity/step.png');

    f = true;
end

function f = temperature_gradient(data, savepath, variable_names, time_indices)
    % Validate inputs

    if length(variable_names) ~= 2
        error('Exactly two variable names must be provided in the variable_names list.');
    end
    
    % Initialize a figure
    figure;
    hold on;
    
    % Loop through each time index
    for t = time_indices

        % Initialize arrays for both variables
        indices1 = [];
        values1 = [];
        indices2 = [];
        values2 = [];
        
        % Extract data for the first variable
        for i = 1:length(data.Properties.VariableNames)
            var_name = data.Properties.VariableNames{i};
            if startsWith(var_name, [variable_names{1}, '_'])
                % Extract the index from the variable name
                index = str2double(regexp(var_name, '\d+', 'match'));
                value = data{t, var_name};
                indices1 = [indices1, index];
                values1 = [values1, value];
            end
        end
        
        % Extract data for the second variable
        for i = 1:length(data.Properties.VariableNames)
            var_name = data.Properties.VariableNames{i};
            if startsWith(var_name, [variable_names{2}, '_'])
                % Extract the index from the variable name
                index = str2double(regexp(var_name, '\d+', 'match'));
                value = data{t, var_name};
                indices2 = [indices2, index];
                values2 = [values2, value];
            end
        end
        
        % Sort indices and corresponding values for both variables
        [indices1, sortIdx1] = sort(indices1);
        values1 = values1(sortIdx1);
        [indices2, sortIdx2] = sort(indices2);
        values2 = values2(sortIdx2);

        % Plot the data for both variables
        plot(indices1+1, values1, '-o', 'DisplayName', [variable_names{1}], 'color', '#333CF5', 'LineWidth', 0.8);
        plot(indices2+1, values2, '-s', 'DisplayName', [variable_names{2}], 'color', '#FFB80E', 'LineWidth', 0.8);
    end

    % Customize the plot
    xlabel('Index');
    ylabel('Temperature [K]');
    legend('T','T_g');
    legend('Location', 'best');
    xlim([1, 37]);
    grid on;
    hold off;
    
    % Save the plot if a savepath is provided
    if savepath
        saveas(gcf, [savepath, 'Temp_comparison.png']);
    end
    
    f = true;
end

function f = plot_propagation(data, var, unit, indices, savepath)
    numElements = length(indices);
    
    fig = figure('Visible', 'on');
    hold on;
    ylabel([var, ' [', unit, ']']);
    xlabel('Time [s]');
    grid on;
 
    % Define line styles and colors to cycle through
    lineStyles = {'-', '--', ':', '-.'};
    colors = lines(numElements); % Use MATLAB's default color palette
    legendEntries = {}; % To store legend labels
    
    for i = 1:numElements
        yData = data.([var, '_', num2str(indices(i))]);
        yDataFiltered = yData(350:600); 
        timeFiltered = linspace(0, length(yDataFiltered), length(yDataFiltered)); 
        
        % Choose line style and color
        lineStyle = lineStyles{mod(i-1, length(lineStyles)) + 1};
        color = colors(i, :);
        
        % Plot the data with specified line style, color, and width
        plot(timeFiltered, yDataFiltered, 'LineStyle', lineStyle, 'Color', color, 'LineWidth', 1.5);
        
        % Add entry for legend
        legendEntries{end+1} = [var, '_', '{',num2str(indices(i)+1),'}']; 
    end

    legend(legendEntries); % Add the legend once after the loop
    hold off;

    if savepath
       saveas(fig, [savepath,'.png']);
    end

    f = true;
end

function f = plot_double_hex(data, data2, var, savepath, end_time)
    numElements = length(var);
    
    i = 1;
    while true
        if i > numElements
            break;
        end

        currentVar = var{i}{1}; % Variable name
        unit = var{i}{2};       % Unit string

        yData = data.(currentVar);
        yData2 = data2.(currentVar);
        
        yDataFiltered = yData(1:end_time); 
        yData2Filtered = yData2(1:end_time);
        time = linspace(0, length(yDataFiltered), length(yDataFiltered));
        time2 = linspace(0, length(yData2Filtered), length(yData2Filtered));
        % correct for 0 indexing in casadi
        parts = split(currentVar, '_');
        segmentVar = parts{1};
        segmentIndex = str2double(parts{2}) + 1;
        segmentIndex = num2str(segmentIndex);

        fig = figure('Visible', 'on');

        hold on
        plot(time, yDataFiltered, '-', 'color', '#333CF5', 'LineWidth',0.8);
        plot(time2, yData2Filtered, '--', 'color', '#FFB80E', 'LineWidth',0.8);
        ylabel([segmentVar, '_{', segmentIndex,'} ', '[', unit, ']']);
        xlabel('Time [s]');
        xlim([0,end_time]);
        legend('Positive step', 'Negative step');
        grid on;
        hold off;

        if savepath
            if segmentVar == 'm'
                saveas(fig, [savepath,'massflow', '_', segmentIndex,'.png']);
            else
                saveas(fig, [savepath, segmentVar, '_', segmentIndex,'.png']);
            end
        end

        %close(fig);
    
        i = i + 1;
    end


    f = true;
end

function f = plot_double(data, data2, var, savepath, end_time)
    numElements = length(var);
    
    i = 1;
    while true
        if i > numElements
            break;
        end

        currentVar = var{i}{1}; % Variable name
        unit = var{i}{2};       % Unit string

        yData = data.(currentVar);
        yData2 = data2.(currentVar);
        
        yDataFiltered = yData(1:end_time); 
        yData2Filtered = yData2(1:end_time);
        time = linspace(0, length(yDataFiltered), length(yDataFiltered));
        time2 = linspace(0, length(yData2Filtered), length(yData2Filtered));
        % correct for 0 indexing in casadi
           
        fig = figure('Visible', 'on');

        hold on
        plot(time, yDataFiltered, '-', 'color', '#333CF5', 'LineWidth',0.8);
        plot(time2, yData2Filtered, '--', 'color', '#FFB80E', 'LineWidth',0.8);
        
        xlabel('Time [s]');
        xlim([0,end_time]);
        legend('Location','best');
        legend('T_g^0 +10K', 'T_g^0 -10K');
        grid on;

        if strcmp(currentVar, 'TU_real')  
           ylabel(['TU_{real}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'm_cold_out')
            ylabel(['m_{cold}^{out}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'T_cold_out')
            ylabel(['T_{cold}^{out}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'm_in')
            ylabel(['m_{p}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hS')
            ylabel(['h_S', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hT')
            ylabel(['h_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'pS')
            ylabel(['p_s', ' [', unit, ']']);
        elseif strcmp(currentVar, 'pT')
            ylabel(['p_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'Ts')
            ylabel(['T_s', ' [', unit, ']']);
            %plot(time, 682*ones(1, end_time), 'LineStyle','--', 'Color','black', 'DisplayName','Setpoint');
        elseif strcmp(currentVar, 'P')
            ylabel(['P', ' [', unit, ']']);
            plot(time, 11500*ones(1, end_time), 'LineStyle','--', 'Color','black', 'DisplayName','Setpoint');
        elseif strcmp(currentVar, 'TT')
            ylabel(['T_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_zm')
            ylabel(['z_m', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_mg')
            ylabel(['m_g^0', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_zV')
            ylabel(['z_v', ' [', unit, ']']);
        elseif strcmp(currentVar, 'HT')
            ylabel(['H_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hU')
            ylabel(['h_U', ' [', unit, ']']);
        elseif strcmp(currentVar, 'MT')
            ylabel(['M_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'Ms')
            ylabel(['M_s', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TBs')
            ylabel(['T_s^{sat}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TBT')
            ylabel(['T_T^{sat}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TBU')
            ylabel(['T_U^{sat}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TU')
            ylabel(['T_U', ' [', unit, ']']);
        elseif strcmp(currentVar, 'rho_s')
            ylabel(['\rho_s', ' [', unit, ']']);
        elseif strcmp(currentVar, 'rho_T')
            ylabel(['\rho_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_m_in_SP')
            ylabel(['m_{p}^{SP}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hc')
            ylabel(['h_c', ' [', unit, ']']);
        else 
           ylabel([currentVar, ' [', unit, ']']);
        end


        hold off;

        if savepath
            if currentVar == 'm'
                saveas(fig, [savepath,'massflow', '.png']);
            elseif currentVar == 'P'
                saveas(fig, [savepath,'Power', '.png']);
            elseif strcmp(currentVar, 'hS')
                saveas(fig, [savepath,'h_S', '.png']);
            elseif strcmp(currentVar, 'hT')
                saveas(fig, [savepath,'h_T', '.png']);
            elseif strcmp(currentVar, 'hc')
                saveas(fig, [savepath,'h_c', '.png']);
            else
                saveas(fig, [savepath, currentVar, '.png']);
            end
        end

        %close(fig);
    
        i = i + 1;
    end


    f = true;
end

function f = plot_triple(data, data2, data3, var, savepath, end_time)
    numElements = length(var);
    
    i = 1;
    while true
        if i > numElements
            break;
        end

        currentVar = var{i}{1}; % Variable name
        unit = var{i}{2};       % Unit string

        yData = data.(currentVar);
        yData2 = data2.(currentVar);
        yData3 = data3.(currentVar);
        
        yDataFiltered = yData(1:end_time); 
        yData2Filtered = yData2(1:end_time);
        yData3Filtered = yData3(1:end_time);

        time = linspace(0, length(yDataFiltered), length(yDataFiltered));
        time2 = linspace(0, length(yData2Filtered), length(yData2Filtered));
        time3 = linspace(0, length(yData3Filtered), length(yData3Filtered));
  
        fig = figure('Visible', 'on');

        hold on
        plot(time, yDataFiltered, '-', 'color', '#333CF5', 'LineWidth',0.8);
        plot(time2, yData2Filtered, '--', 'color', '#FFB80E', 'LineWidth',0.8);
        plot(time3, yData3Filtered, "-.", 'color', '#BC19BF', 'LineWidth',0.8);

        if strcmp(currentVar, 'TU_real')  
           ylabel(['TU_{real}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'm_cold_out')
            ylabel(['m_{cold}^{out}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'T_cold_out')
            ylabel(['T_{cold}^{out}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'm_in')
            ylabel(['m_{p}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hS')
            ylabel(['h_S', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hT')
            ylabel(['h_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'pS')
                ylabel(['p_s', ' [', unit, ']']);
        elseif strcmp(currentVar, 'pT')
                ylabel(['p_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'Ts')
                ylabel(['T_s', ' [', unit, ']']);
                plot(time, 682*ones(1, end_time), 'LineStyle','--', 'Color','black', 'DisplayName','Setpoint','LineWidth',0.8);
        elseif strcmp(currentVar, 'TT')
                ylabel(['T_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'dummy_zm')
                ylabel(['z_m', ' [', unit, ']']);
        elseif strcmp(currentVar, 'HT')
                ylabel(['H_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'Hs')
                ylabel(['H_s', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hU')
                ylabel(['h_U', ' [', unit, ']']);
        elseif strcmp(currentVar, 'MT')
                ylabel(['M_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'Ms')
                ylabel(['M_s', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TBs')
                ylabel(['T_s^{sat}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TBT')
                ylabel(['T_T^{sat}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TBU')
                ylabel(['T_U^{sat}', ' [', unit, ']']);
        elseif strcmp(currentVar, 'TU')
                ylabel(['T_U', ' [', unit, ']']);
        elseif strcmp(currentVar, 'hc')
            ylabel(['h_c', ' [', unit, ']']);
        elseif strcmp(currentVar, 'rho_s')
            ylabel(['\rho_s', ' [', unit, ']']);
        elseif strcmp(currentVar, 'rho_T')
            ylabel(['\rho_T', ' [', unit, ']']);
        elseif strcmp(currentVar, 'Qc')
            ylabel(['Q_c', ' [', unit, ']']);  
        elseif strcmp(currentVar, 'dummy_m_in_SP')
            ylabel(['m_{p}^{SP}', ' [', unit, ']']);
            
        else 
           ylabel([currentVar, ' [', unit, ']']);
        end

        xlabel('Time [s]');
        xlim([0,end_time]);
        legend('Location','best');
        legend('Feedback + FF', 'Feedforward', 'Feedback');
        grid on;
        hold off;

        if savepath
            if currentVar == 'm'
                saveas(fig, [savepath,'massflow', '.png']);
            elseif currentVar == 'P'
                saveas(fig, [savepath,'Power', '.png']);
            elseif strcmp(currentVar, 'hS')
                saveas(fig, [savepath,'h_S', '.png']);
            elseif strcmp(currentVar, 'hT')
                saveas(fig, [savepath,'h_T', '.png']);
            elseif strcmp(currentVar, 'hc')
                saveas(fig, [savepath,'h_c', '.png']);
            else
                saveas(fig, [savepath, currentVar, '.png']);
            end
        end

    
        i = i + 1;
    end


    f = true;
end

function f = plot_two(data, var, savepath, end_time)

    currentVar1 = var{1}; % First variable name
    currentVar2 = var{2}; % Second variable name

    yData = data.(currentVar1);
    yData2 = data.(currentVar2);

    yDataFiltered = yData(1:end_time); 
    yData2Filtered = yData2(1:end_time);

    time = linspace(0, length(yDataFiltered), length(yDataFiltered));
    time2 = linspace(0, length(yData2Filtered), length(yData2Filtered));

    fig = figure('Visible', 'on');

    %% First subplot
    subplot(2,1,1);
    plot(time, yDataFiltered, '-', 'color', '#333CF5', 'LineWidth', 0.8);
    hold on;

    % Crossings at 0
    crossZero1 = find(yDataFiltered(1:end-1).*yDataFiltered(2:end) < 0);
    % Crossings at 1 (value goes above or below 1)
    crossOne1 = find((yDataFiltered(1:end-1) - 1).*(yDataFiltered(2:end) - 1) < 0);

    % Plot markers
    plot(time(crossZero1), yDataFiltered(crossZero1), 'ko', 'MarkerSize', 5, 'DisplayName', 'Cross 0');
    plot(time(crossOne1), yDataFiltered(crossOne1), 'ro', 'MarkerSize', 5, 'DisplayName', 'Cross 1');

    ylabel([currentVar1, ' [unit]']);
    xlabel('Time [s]');
    xlim([0, end_time]);
    title(['Time Series of ', currentVar1]);
    legend('Location','best');
    grid on;
    hold off;

    %% Second subplot
    subplot(2,1,2);
    plot(time2, yData2Filtered, '--', 'color', '#FFB80E', 'LineWidth', 0.8);
    hold on;

    % Crossings at 0
    crossZero2 = find(yData2Filtered(1:end-1).*yData2Filtered(2:end) < 0);
    % Crossings at 1
    crossOne2 = find((yData2Filtered(1:end-1) - 1).*(yData2Filtered(2:end) - 1) < 0);

    % Plot markers
    plot(time2(crossZero2), yData2Filtered(crossZero2), 'ko', 'MarkerSize', 5, 'DisplayName', 'Cross 0');
    plot(time2(crossOne2), yData2Filtered(crossOne2), 'ro', 'MarkerSize', 5, 'DisplayName', 'Cross 1');

    ylabel([currentVar2, ' [unit]']);
    xlabel('Time [s]');
    xlim([0, end_time]);
    title(['Time Series of ', currentVar2]);
    legend('Location','best');
    grid on;
    hold off;

    %% Save if path provided
    if savepath
        saveas(fig, [savepath, currentVar1, '_and_', currentVar2, '.png']);
    end

    f = true;
end

function f = plot_cross_beta(data, var, savepath, end_time)

    currentVar1 = var{1}; % First variable name

    yData = data.(currentVar1);
    yDataFiltered = yData(1:end_time); 
    time = linspace(0, length(yDataFiltered), length(yDataFiltered));

    fig = figure('Visible', 'on');

    %% First plot
    p1 = plot(time, yDataFiltered, '-', 'color', '#333CF5', 'LineWidth', 0.8);
    hold on;

    % Hide first plot from legend
    p1.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Crossings at 0
    crossZero1 = find(yDataFiltered(1:end-1).*yDataFiltered(2:end) < 0);
    % Crossings at 1
    crossOne1 = find((yDataFiltered(1:end-1) - 1).*(yDataFiltered(2:end) - 1) < 0);

    % Plot crossing markers
    plot(time(crossZero1), yDataFiltered(crossZero1), 'ko', 'MarkerSize', 5, 'DisplayName', 'Switch \beta = 0');
    plot(time(crossOne1), yDataFiltered(crossOne1), 'ro', 'MarkerSize', 5, 'DisplayName', 'Switch \beta = 1');

    ylabel('\beta_{28} [-]');
    xlabel('Time [s]');
    xlim([0, end_time]);
    legend('Location','southeast');
    grid on;
    hold off;

    %% Save if path provided
    if savepath
        saveas(fig, [savepath, currentVar1, '_cross', '.png']);
    end

    f = true;
end

function f = plot_cross(data, var, savepath, end_time, crossVars, crossingTargets)

    currentVar1 = var{1}; % First variable name
    unit = var{2};        % Unit string

    yData = data.(currentVar1);
    yDataFiltered = yData(1:end_time); 
    time = linspace(0, length(yDataFiltered), length(yDataFiltered));

    fig = figure('Visible', 'on');

    %% First plot
    p1 = plot(time, yDataFiltered, '-', 'color', '#333CF5', 'LineWidth', 0.8);
    hold on;
    % Hide first plot from legend
    p1.Annotation.LegendInformation.IconDisplayStyle = 'off';

    %% Initialize arrays for crossings
    allCrossings0 = []; % Crossings at 0
    allCrossings1 = []; % Crossings at 1

    %% Loop over crossing variables
    for i = 1:length(crossVars)
        crossVarName = crossVars{i};
        targetValue = crossingTargets(i);

        yCross = data.(crossVarName);
        yCrossFiltered = yCross(1:end_time);

        crossings = find((yCrossFiltered(1:end-1) - targetValue) .* (yCrossFiltered(2:end) - targetValue) < 0);

        % Save crossings according to target
        if targetValue == 0
            allCrossings0 = [allCrossings0; crossings(:)];
        elseif targetValue == 1
            allCrossings1 = [allCrossings1; crossings(:)];
        end
    end

    %% Plot all 0-crossings
    if ~isempty(allCrossings0)
        p0 = plot(time(allCrossings0), yDataFiltered(allCrossings0), 'ko', 'MarkerSize', 4, 'LineWidth', 1.2,'DisplayName', 'Switch \beta = 0');
    end

    %% Plot all 1-crossings
    if ~isempty(allCrossings1)
        p1 = plot(time(allCrossings1), yDataFiltered(allCrossings1), 'ro', 'MarkerSize', 4, 'LineWidth', 1.2,'DisplayName', 'Switch \beta = 1');
    end

    ylabel(['M_{11}', ' [', unit, ']']);
    xlabel('Time [s]');
    xlim([0, end_time]);
    legend('Location', 'best');
    grid on;
    hold off;

    %% Save if path provided
    if savepath
        saveas(fig, [savepath, currentVar1, '_cross', '.png']);
    end

    f = true;
end

function rmse = rmse_deviation(data, var, end_time, setpoint)

    variableName = var{1}; % Extract variable name
    yData = data.(variableName);
    yDataFiltered = yData(1:end_time);

    error = yDataFiltered - setpoint;
    rmse = sqrt(mean(error.^2));
end




