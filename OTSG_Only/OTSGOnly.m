function [x_keys, z_keys, sol] = OTSGOnly(n, x_0, z_0, constants)
    import casadi.*

    %% constants
    R       = constants(1);  % [m3*bar/K/mol]
    cpS     = constants(2);  % [kJ/kg/K] 
    cpW     = constants(3);  % [kJ/kg/K] 
    cpG     = constants(4);  % [kJ/kg/K] 
    T_Ref   = constants(5);  % [K] 
    TB_Ref  = constants(6);  % [K] 
    T0g     = constants(7);  % [K] 
    Tp      = constants(8);  % [K]
    mG      = constants(9);  % [kg/s] 
    dHvap0  = constants(10); % [kJ/kg] 
    Us      = constants(11); % [kW/m2/K]
    Ue      = constants(12); % [kW/m2/K]
    Mw      = constants(13); % [kg/mol] 
    Cvd     = constants(14); % [kg/bar]
    p_in    = constants(15); % [bar]
    p_out   = constants(16); % [bar]
    V_tot   = constants(17); % [m3]
    rho_Ref = constants(18); % [kg/m3]
    p_Ref   = constants(19); % [bar]
    k_p     = constants(20); % [1/bar]
    A_anto  = constants(21); % [-]
    B_anto  = constants(22); % [K]
    C_anto  = constants(23); % [K]
    A_tot   = constants(24); % [m2]
    A_corr  = constants(25);
    
    % Segment dependent constants
    V     = V_tot/n;
    A     = A_tot/n; 


    %% Solver initialization
    % Define variables
    m_in   = SX.sym('m_in');
    m      = SX.sym('m', n);
    T      = SX.sym('T', n);
    h      = SX.sym('h', n);
    M      = SX.sym('M', n);
    H      = SX.sym('H', n);
    Tg     = SX.sym('Tg', n);
    mg     = SX.sym('mg', n);
    Q      = SX.sym('Q', n);
    p      = SX.sym('p', n);
    rho    = SX.sym('rho', n);
    beta   = SX.sym('beta', n);
    TB     = SX.sym('TB', n);
  
    % Set states  x: differential, z: algebraic
    x = [];            
    for k=1:n
        x = [x; M(k);H(k)];
    end

    z = [m_in];
    for k=1:n
        z = [z; m(k);T(k);h(k);Tg(k);mg(k);Q(k);p(k);rho(k);beta(k);TB(k)];
    end

    % Initialize equatuions
    Alg  = [];
    diff = [];

    % Set function output variables
    x_keys = x;    
    z_keys = z;

    %% Additional inlet equations
    h_in  = cpW*(Tp-T_Ref);        % Note: Assuming input is liquid state
    init1 = m_in - Cvd*(p_in-p(1));
    Alg   = [Alg;init1];

    %% i = 1
    dHvap = dHvap0 + (cpW-cpS)*(TB_Ref-TB(1)); 
    
    % Switch logic
    cond1 = beta(1) >= 1;
    cond2 = beta(1) <= 0;

    % 1: Steam, 2: Two-phase, 3: Liquid
    Enth1 = cpW*(TB(1)-T_Ref) + dHvap + cpS*(T(1)-TB(1)) - h(1); 
    Enth2 = T(1) - TB(1);                                          
    Enth3 = cpW*(T(1)-T_Ref) - h(1);                              

    % Compute U using numeric A_corr
    U = U_Calc(beta(1), Tg(1), A_corr); % Pass numeric A_corr
    UA = A*U;
    
    VL    = (1-beta(1))*M(1)/rho(1);
    Pres1 = M(1)*R*T(1)/(Mw*V) - p(1);
    Pres2 = p(1)*(V-VL) - beta(1)*M(1)*R*T(1)/(Mw);
    Pres3 = rho(1) - M(1)/V;

    % Equations
    dMdt = m_in - m(1);
    dHdt = m_in*h_in - m(1)*h(1) + Q(1);

    alg1  = mg(1)*cpG*(Tg(2)-Tg(1)) - Q(1);
    alg2  = UA*((Tg(2)-T(1))/2+(Tg(1)-Tp)/2) - Q(1);
    alg3  = m(1) - Cvd*(p(1)-p(2)); 
    alg4  = if_else(cond1, Enth1, if_else(cond2, Enth3, Enth2));
    alg5  = M(1)*h(1) - H(1);
    alg6  = mg(1) - mg(2);
    alg7  = if_else(cond1, Pres1, if_else(cond2, Pres3, Pres2));
    alg8  = p(1) - 1/(k_p*rho_Ref)*(rho(1)-rho_Ref) - p_Ref;
    alg9  = cpW*(TB(1)-T_Ref) + beta(1)*dHvap - h(1);
    alg10 = 10^(A_anto-(B_anto/(TB(1)+C_anto))) - p(1);

    Alg  = [Alg;alg1;alg2;alg3;alg4;alg5;alg6;alg7;alg8;alg9;alg10];
    diff = [diff;dMdt;dHdt];

    %% i = 2 to n-1
    for k=2:n-1
        dHvap = dHvap0 + (cpW-cpS)*(TB_Ref-TB(k)); 
        % Switch logic
        cond1 = beta(k) >= 1;
        cond2 = beta(k) <= 0;
    
        Enth1 = cpW*(TB(k)-T_Ref) + dHvap + cpS*(T(k)-TB(k)) - h(k); 
        Enth2 = T(k) - TB(k);                                          
        Enth3 = cpW*(T(k)-T_Ref) - h(k);                             

        U = U_Calc(beta(k), Tg(k), A_corr); 
        UA = A*U;

        VL    = (1-beta(k))*M(k)/rho(k);
        Pres1 = M(k)*R*T(k)/(Mw*V) - p(k);
        Pres2 = p(k)*(V-VL) - beta(k)*M(k)*R*T(k)/(Mw);
        Pres3 = rho(k) - M(k)/V;
 
        % Equations
        dMdt = m(k-1)-m(k);
        dHdt = m(k-1)*h(k-1) - m(k)*h(k) + Q(k);

        alg1  = mg(k)*cpG*(Tg(k+1)-Tg(k)) - Q(k);
        alg2  = UA*((Tg(k+1)-T(k))/2+(Tg(k)-T(k-1))/2) - Q(k);
        alg3  = m(k) - Cvd*(p(k)-p(k+1));
        alg4  = if_else(cond1, Enth1, if_else(cond2, Enth3, Enth2));
        alg5  = M(k)*h(k) - H(k);
        alg6  = mg(k) - mg(k+1);
        alg7  = if_else(cond1, Pres1, if_else(cond2, Pres3, Pres2));
        alg8  = p(k) - 1/(k_p*rho_Ref)*(rho(k)-rho_Ref) - p_Ref;
        alg9  = cpW*(TB(k)-T_Ref) + beta(k)*dHvap - h(k);
        alg10 = 10^(A_anto-(B_anto/(TB(k)+C_anto))) - p(k);

        Alg  = [Alg;alg1;alg2;alg3;alg4;alg5;alg6;alg7;alg8;alg9;alg10];
        diff = [diff;dMdt;dHdt];
    end

    %% i = n
    dHvap = dHvap0 + (cpW-cpS)*(TB_Ref-TB(n)); 
    % Switch logic
    cond1 = beta(n) >= 1;
    cond2 = beta(n) <= 0;

    Enth1 = cpW*(TB(n)-T_Ref) + dHvap + cpS*(T(n)-TB(n)) - h(n); 
    Enth2 = T(n) - TB(n);                                          
    Enth3 = cpW*(T(n)-T_Ref) - h(n);                             

    U = U_Calc(beta(n), Tg(n), A_corr); 
    UA = A*U;

    VL    = (1-beta(n))*M(n)/rho(n);
    Pres1 = M(n)*R*T(n)/(Mw*V) - p(n);
    Pres2 = p(n)*(V-VL) - beta(n)*M(n)*R*T(n)/(Mw);
    Pres3 = rho(n) - M(n)/V;

    % Equations
    dMdt = m(n-1)-m(n);
    dHdt = m(n-1)*h(n-1) - m(n)*h(n) + Q(n);

    alg1  = mg(n)*cpG*(T0g-Tg(n)) - Q(n);
    alg2  = UA*((T0g-T(n))/2+(Tg(n)-T(n-1))/2) - Q(n);
    alg3  = m(n) - Cvd*(p(n)-p_out);
    alg4  = if_else(cond1, Enth1, if_else(cond2, Enth3, Enth2));
    alg5  = M(n)*h(n) - H(n);
    alg6  = mg(n) - mG;
    alg7  = if_else(cond1, Pres1, if_else(cond2, Pres3, Pres2));
    alg8  = p(n) - 1/(k_p*rho_Ref)*(rho(n)-rho_Ref) - p_Ref;
    alg9  = cpW*(TB(n)-T_Ref) + beta(n)*dHvap - h(n);
    alg10 = 10^(A_anto-(B_anto/(TB(n)+C_anto))) - p(n);

    Alg  = [Alg;alg1;alg2;alg3;alg4;alg5;alg6;alg7;alg8;alg9;alg10];
    diff = [diff;dMdt;dHdt];
    
    %% Solver
    dae     = struct;
    dae.x   = x;     % Differential states
    dae.z   = z;     % Algebraic states
    dae.ode = diff;  % Differential equations
    dae.alg = Alg;   % Algebraic equations

    opts = struct('tf', 1, 'abstol', 1e-9, 'reltol', 1e-9, ...
        'max_num_steps', 10000);
    
    F = integrator('F', 'idas', dae, opts);
   
    sol = F('x0', x_0, 'z0', z_0);
end