clear
clc
addpath('casadi-3.6.6')


save = {'test1.csv', 'test2.csv'};

T0g = 443.338+273.15;

Tg_step_N = [(T0g)*ones(1,200), (T0g-10)*ones(1,3800)];
Tg_step_P = [(T0g)*ones(1,200), (T0g+10)*ones(1,3800)];

applied_step = {Tg_step_N, Tg_step_P};
step_index  = {7, 7};

for y=1:2
% Set filepath for storage, use false to discard results
filepath = save{y};
load_from_file = 'load_file.csv';
Hex_only = false;
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
Mw      = 18*10^(-3);       % [kg/mol] molar weight water
p_pump  = 29;
p_in    = 24;               % [bar] Cold side inlet pressure 
pS      = 23;               % [bar] Cold side outlet pressure
pT      = 22;               % [bar] Turbine inlet pressure
V_tot   = 3.92;             % [m3] OTSG volume
A_tot   = 739.4;            % [m2] OTSG Area
A_corr  = 5389.79/A_tot;    % [-] Fin correction on hot side
rho_Ref = 1000;             % [kg/m3] Reference density water
p_Ref   = 1;                % [bar] Reference pressure water
k_p     = 4.58*10^(-5);     % [1/bar] Compresibility factor water
A_anto  = 3.55959;          % [-] Antoine Coeff (379. to 573K)
B_anto  = 643.748;          % [K] Antoine Coeff (379. to 573K)
C_anto  = -198.043;         % [K] Antoine Coeff (379. to 573K)
A_anto_c= 4.6543;           % [-] Antoine Coeff (255.9 to 373K)
B_anto_c= 1435.264;         % [K] Antoine Coeff (255.9 to 373K)
C_anto_c= -64.848;          % [K] Antoine Coeff (255.9 to 373K)
Cvd     = 21.9/(2*(p_in-pS))*(n+1); % [kg/bar] 
Kv      = 21.9/(2*0.9*(pS-pT));     % [kg/bar] valve coeff for valve after OTSG
zV      = 0.9;              % valve opening for valve after OTSG
Vs      = 0.5;              % [m3]  Steam buffer volume
phiD    = 13;               % [kgK^1/2/(s bar)] Turbine massflow coefficient
Mg      = 2.5022e+03;       % Mg0*P0/w0; 
Dg      = 661.9580;         % Dg0*P0/w0; 
w0      = 50;               % [Hz] Reference grid frequency
L       = 11;               % [kW] denoted as PL in code 
p_c       = 0.0358;         % [bar] outlet pressure boundary
eta       = 0.9;            % [-] Turbine efficency factor
VT        = 0.5; 

p_guess  = linspace(p_in, pS, n+2);
p_guess  = p_guess(2:n+1);    % [bar]
k_m      = 21.9/(2*0.5*(p_pump-p_guess(1))); % [kg/bar]

constants = [R, cpS, cpW, cpG, T_Ref, TB_Ref, T0g, Tp, mG, dHvap0, ...
    Mw, Cvd, V_tot, rho_Ref, p_Ref, k_p, A_anto, B_anto, C_anto,...
    A_tot, A_corr, Kv, zV, Vs, phiD, Mg, Dg, w0, L, ...
    A_anto_c, B_anto_c, C_anto_c, eta, VT, k_m, p_pump, p_c];

%% Set guesses 
mP         = 21.9/2;          % [kg/s]
beta_guess = -0.4;            % [-]
M_guess    = V_tot*1000/(n);  % [kg]
T_guess    = Tp;              % [K]
h_guess    = cpW*(Tp-T_Ref);  % [kJ/kg]
H_guess    = M_guess*h_guess; % [kJ]
m_guess    = mP;              % [kg/s]
Tg_guess   = T0g;             % [K]
mg_guess   = mG;              % [kg/s]
U_guess    = 0.4;             % [kW/K]
Q_guess    = U_guess*A_tot/n*(Tg_guess-T_guess);  % [kW]
rho_guess  = 1000;            % [kg/m3]


TB_guess = TB_Ref;            % [K]

m_s_guess   = mP;
Ts_guess    = 683;            % [K]
TBs_guess   = TB_Ref;         % [K]
hs_guess    = cpW*(TBs_guess-T_Ref) + dHvap0 + cpS*(Ts_guess-TBs_guess); % [kJ/kg]
beta_s_guess= (hs_guess-cpW*(TBs_guess-T_Ref))/dHvap0; % [-]
rho_s_guess = rho_guess;      % [kg/m3]
Ms_guess    = pS*Vs*Mw/(R*Ts_guess);  % [Kg]
Hs_guess    = Ms_guess*hs_guess;      % [kJ]

m_T_guess   = mP;             % [kg/s]
TT_guess    = Ts_guess;       % [K]
TBT_guess   = TB_Ref;
hT_guess    = cpW*(TBT_guess-T_Ref) + dHvap0 + cpS*(TT_guess-TBT_guess); % [kJ/kg]
beta_T_guess= (hT_guess-cpW*(TBT_guess-T_Ref))/dHvap0; % [-]
rho_T_guess = rho_guess;
MT_guess    = pT*VT*Mw/(R*TT_guess);  % [Kg]
HT_guess    = MT_guess*hT_guess;      % [kJ]

W_guess     = -16500/2;       % [kW]
P_guess     = -W_guess;       % [kW]
w_guess     = 50;             % [Hz]

m_U_guess   = mP;
TU_guess    = 240;
TU_real_guess= 273.15+27;                              % [K]
TBU_guess   = TU_real_guess;                           % [K]
hU_guess    = cpW*(TBU_guess-T_Ref) + dHvap0 + cpS*(TU_guess-TBU_guess); % [kJ/kg]
beta_U_guess= 0.9;



Tc_guess    = TU_real_guess;                               % [K]
m_c_guess   = mP;                                      % [kg/s]
pc_guess    = 0.0358;                                  % [bar]
hc_guess    = cpW*(Tc_guess-T_Ref);                    % [kJ/kg]
Qc_guess    = m_U_guess*hU_guess - m_c_guess*hc_guess; % [kW]
Mb_guess    = 10000;                                   % [kg] buffer tank 

% Control
control = [];
%Flow controller
m_in_SP_0   = mP;
z_m_0     = 0.5;
KI_m_in   = 0.009699;
Ts_SP     = 682; % CHANGE THIS?
Kc_Ts     = -0.03036;
tauI_Ts   = 236;

control = [control; m_in_SP_0; z_m_0; KI_m_in; Ts_SP; Kc_Ts; tauI_Ts];

dummy_zm_guess = z_m_0;
dummy_m_in_SP  = m_in_SP_0;
interrmIn_guess = 0;
interrTs_guess = 0;


% Set order of guesses to match order in solver
if ~isempty(load_from_file) && Hex_only
    loadTable = readtable(load_from_file);
    lastRow = loadTable{end, :};  % Want the last instance in time
    
    % OTSG variables
    num_x = 2*n;    % Number of x variables
    num_z = 10*n+1; % Number of z variables

    % Reorder: x and z are swapped in key_names
    reordered_data = [lastRow(num_z+1:end), lastRow(1:num_z)]; 

    % Assign to initial conditions
    x_0 = reordered_data(1:num_x); 
    z_0 = reordered_data(num_x+1:end)'; 
    
    % Variables after OTSG 
    x_0 = [x_0, Ms_guess, Hs_guess, MT_guess, HT_guess, w_guess, ...
        Mb_guess];

    z_0 = [z_0; m_s_guess; Ts_guess; TBs_guess; pS; hs_guess; beta_s_guess; ...
        rho_s_guess; m_T_guess; TT_guess; TBT_guess; pT; hT_guess; ...
        beta_T_guess; rho_T_guess; W_guess; P_guess; m_U_guess; TU_guess;...
        TU_real_guess; TBU_guess; hU_guess; beta_U_guess; ...
        Tc_guess; m_c_guess; hc_guess; Qc_guess];

    z_0 = [z_0; dummy_zm_guess; dummy_m_in_SP];
    x_0 = [x_0, interrmIn_guess, interrTs_guess];

elseif ~isempty(load_from_file) && ~Hex_only
    loadTable = readtable(load_from_file);
    lastRow = loadTable{end, :};

    % System variables
    num_x = 2*n+8;     % Number of x variables
    num_z = 10*n+1+28; % Number of z variables (hex + inlet + other)

    reordered_data = [lastRow(num_z+1:end), lastRow(1:num_z)]; 

    % Assign to initial conditions
    x_0 = reordered_data(1:num_x); 
    z_0 = reordered_data(num_x+1:end)'; 

else
    x_0 = repmat([M_guess; H_guess], n, 1);
    z_0 = [m_guess]; 
    for k=1:n
        z_0 = [z_0; m_guess; T_guess; h_guess; Tg_guess; mg_guess; Q_guess; ...
            p_guess(k); rho_guess; beta_guess;TB_guess];
    end
    z_0 = [z_0; m_s_guess; Ts_guess; TBs_guess; pS; hs_guess; beta_s_guess; ...
        rho_s_guess; m_T_guess; TT_guess; TBT_guess; pT; hT_guess; ...
        beta_T_guess; rho_T_guess; W_guess; P_guess; m_U_guess; TU_guess...
        TU_real_guess; TBU_guess; hU_guess; beta_U_guess; ...
        Tc_guess; m_c_guess; hc_guess; Qc_guess];

    x_0 = [x_0, Ms_guess, Hs_guess, MT_guess, HT_guess, w_guess, ...
         Mb_guess];

    z_0 = [z_0; dummy_zm_guess; dummy_m_in_SP];
    x_0 = [x_0, interrmIn_guess; interrTs_guess];
end


%% Run
T = 4000; % [s] Total simulation time
N = T;   % [-] Number of steps 

% Inital conditions
x = x_0;
z = z_0;

% Make result arrays
x_save = zeros(N+1, length(x_0)); 
z_save = zeros(N+1, length(z_0)); 
x_save(1, :) = x_0;
z_save(1, :) = z_0;

% Steps

time_save = [0];
running_time = 0;
for i = 1:N

    % Apply steps
    constants(step_index{y}) = applied_step{y}(i);  % Assign corresponding value at index

    % Evaluate system 
    [xf_keys, zf_keys, result] = C3Feedback(n, x, z, constants, control); 
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
end
