clear
clc
addpath('casadi-3.6.6')

% Set filepath for storage, use false to discard results
filepath = 'test1.csv';
load_from_file = false;

% Set number of segments
n = 37; 

%% Constants
R       = 8.314462618*1e-5; % [m3*bar/K/mol] gas constant 
cpS     = 2.43;             % [kJ/kg/K] heat capasity steam 
cpW     = 4.24;             % [kJ/kg/K] heat capasity water 
cpG     = 1.02;             % [kJ/kg/K] heat capasity gas   
T_Ref   = 0;                % [K] Reference temperature 
TB_Ref  = 303+273.15;       % [K] Reference boiling point temperature
T0g     = 443.338+273.15;   % [K] Hot side inlet temperature
Tp      = 26.6769+273.15;   % [K] Cold side inlet temperature
mG      = 225.5/2;          % [kg/s] Hot side mass flow
dHvap0  = 1382;             % [kJ/kg] Reference energy of vaporization
Us      = 0.4;              % [kW/m2/K] Heat transfer coeff gas    
Ue      = 0.4;              % [kW/m2/K] Heat transfer coeff liquid 
Mw      = 18*10^(-3);       % [kg/mol] molar weight water
Cvd     = 21.9/(2*1)*(n+1); % [kg/bar] 
p_in    = 24;               % [bar] Cold side inlet pressure
p_out   = 23;               % [bar] Cold side outlet pressure
V_tot   = 3.92;             % [m3] OTSG volume
A_tot   = 739.4;            % [m2] OTSG Area
A_corr  = 5389.79/A_tot;    % Fin correction on hot side
rho_Ref = 1000;             % [kg/m3] Reference density water
p_Ref   = 1;                % [bar] Reference pressure water
k_p     = 4.58*10^(-4);     % [1/bar] Compresibility factor water
A_anto  = 3.55959;          % [-] Antoine Coeff
B_anto  = 643.748;          % [K] Antoine Coeff
C_anto  = -198.043;         % [K] Antoine Coeff

constants = [R, cpS, cpW, cpG, T_Ref, TB_Ref, T0g, Tp, mG, dHvap0, Us, Ue, ...
    Mw, Cvd, p_in, p_out, V_tot, rho_Ref, p_Ref, k_p, A_anto, B_anto, C_anto, A_tot, A_corr];


%% Set guesses 
mP         = 21.9/2;            % [kg/s]
beta_guess = -0.4;            % [-]
M_guess    = V_tot*1000/(n);  % [kg]
T_guess    = Tp;              % [K]
h_guess    = cpW*(Tp-T_Ref);  % [kJ/kg]
H_guess    = M_guess*h_guess; % [kJ]
m_guess    = mP;              % [kg/s]
Tg_guess   = T0g;             % [K]
mg_guess   = mG;              % [kg/s]
Q_guess    = Ue*A_tot/n*(Tg_guess-T_guess);           % [kW]
rho_guess  = 1000;            % [kg/m3]

p_guess  = linspace(p_in, p_out, n+2);
p_guess  = p_guess(2:n+1);    % [bar]
TB_guess = TB_Ref;            % [K]


% Set order of guesses to match order in solver
if load_from_file
    loadTable = readtable(load_from_file);
    lastRow = loadTable{end, :};
    % Reorder data to swap x and z positions
    num_z = 10*n+1; % Number of z variables
    num_x = 2*n; % Number of x variables
    
    % Reorder: x and z are swapped in key_names
    reordered_data = [lastRow(num_z+1:end), lastRow(1:num_z)]; 

    % Assign to initial conditions
    x_0 = reordered_data(1:num_x); 
    z_0 = reordered_data(num_x+1:end); 
else
    x_0 = repmat([M_guess; H_guess], n, 1);
    z_0 = [m_guess]; 
    for k=1:n
        z_0 = [z_0; m_guess; T_guess; h_guess; Tg_guess; mg_guess; Q_guess; ...
            p_guess(k); rho_guess; beta_guess;TB_guess];
    end
end

%% Run
T = 1000; % [s] Total simulation time
N = T;   % [-] Number of steps 

% Inital conditions
x = x_0;
z = z_0;

% Make result arrays
x_save = zeros(N+1, length(x_0)); 
z_save = zeros(N+1, length(z_0)); 
x_save(1, :) = x_0;
z_save(1, :) = z_0;

time_save = [0];
running_time = 0;
for i = 1:N
    % Evaluate system 
    [xf_keys, zf_keys, result] = OTSGOnly(n, x, z, constants); 
    x = result.xf;
    z = result.zf;
    
    % Save state
    x_save(i+1, :) = full(x);
    z_save(i+1, :) = full(z);
    running_time   = running_time + T/N;
    disp(running_time);
    time_save = [time_save, running_time];
end

%% Make result arrays

% Convert keys to strings
zf_keys = cellfun(@(var) char(var.name()), num2cell(zf_keys), 'UniformOutput', false);
xf_keys = cellfun(@(var) char(var.name()), num2cell(xf_keys), 'UniformOutput', false);

% Initialize cell arrays for results
key_names = [zf_keys(:); xf_keys(:)]; % Convert both to column vectors and concatenate
values = [z_save, x_save]; % Combine z_save and x_save in order

% For terminal display
for i = 1:length(key_names)
    data = full(values(:, i));
    data = data(end);
    disp([key_names{i}, ': ', mat2str(round(data,4))])
end

%% Store as CSV

if filepath
    % Create a cell array for storing data
    numKeys = length(key_names);         % Total number of keys
    maxRows = size(values, 1);          % Assuming same number of rows 
    tableData = cell(maxRows + 1, numKeys); % +1 for the header (keys)
    
    % Store the keys as headers in the first row
    for i = 1:numKeys
        tableData{1, i} = key_names{i}; 
    end
    
    % Store the values under each corresponding header
    for i = 1:numKeys
        dataArray = values(:, i); 
        for j = 1:maxRows
            % Store values starting from row 2
            tableData{j + 1, i} = dataArray(j); 
        end
    end
    
    % Convert the cell array to a table
    csvTable = cell2table(tableData(2:end, :), 'VariableNames', tableData(1, :));
    
    % Write the table to a CSV file
    writetable(csvTable, filepath);
end
