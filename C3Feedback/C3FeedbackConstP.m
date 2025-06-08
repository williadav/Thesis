function [x_keys, z_keys, sol] = ...
    C3FeedbackConstP(n, x_0, z_0, constants, control)
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
    MW      = constants(11); % [kg/mol] 
    Cvd     = constants(12); % [kg/bar]
    V_tot   = constants(13); % [m3]
    rho_Ref = constants(14); % [kg/m3]
    p_Ref   = constants(15); % [bar]
    k_p     = constants(16); % [1/bar]
    A_anto  = constants(17); % [-]
    B_anto  = constants(18); % [K]
    C_anto  = constants(19); % [K]
    A_tot   = constants(20); % [m2]
    A_corr  = constants(21); % [-]
    Kv      = constants(22); % [kg/bar] 
    Vs      = constants(23); % [m3]
    phiD    = constants(24);
    Mg      = constants(25);
    Dg      = constants(26);
    w0      = constants(27);
    L       = constants(28);
    A_anto_c= constants(29);
    B_anto_c= constants(30);
    C_anto_c= constants(31);
    eta     = constants(32);
    VT      = constants(33);
    k_m     = constants(34);
    p_pump  = constants(35);
    pc     = constants(36);

    % Controller constants
    m_in_SP_0 = control(1);
    z_m_0     = control(2);
    KI_m_in   = control(3);
    Ts_SP     = control(4);
    Kc_Ts     = control(5);
    tauI_Ts   = control(6);
    pS_SP     = control(7);
    zV_0      = control(8);
    Kc_pS     = control(9);
    tauI_pS   = control(10);
    tauT_pS   = control(11);

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

    Ms     = SX.sym('Ms');
    Hs     = SX.sym('Hs');
    m_s    = SX.sym('m_s');
    Ts     = SX.sym('Ts');
    TBs    = SX.sym('TBs');
    pS     = SX.sym('pS');
    hS     = SX.sym('hS');
    beta_s = SX.sym('beta_s');
    rho_s  = SX.sym('rho_s');

    MT     = SX.sym('MT');
    HT     = SX.sym('HT');
    m_T    = SX.sym('m_T');
    TT     = SX.sym('TT');
    TBT    = SX.sym('TBT');
    pT     = SX.sym('pT');
    hT     = SX.sym('hT');
    beta_T = SX.sym('beta_T');
    rho_T  = SX.sym('rho_T');

    w      = SX.sym('w');
    W      = SX.sym('W');
    P      = SX.sym('P');
    
    m_U    = SX.sym('m_U');
    TU     = SX.sym('TU');
    TU_real= SX.sym('TU_real');
    TBU    = SX.sym('TBU');
    hU     = SX.sym('hU');
    beta_U = SX.sym('beta_U');

    Tc     = SX.sym('Tc');
    m_c    = SX.sym('m_c');
    Qc     = SX.sym('Qc');
    hc     = SX.sym('hc');

    Mb     = SX.sym('Mb');

    % Control
    dummy_zm  = SX.sym('dummy_zm');
    interrmIn = SX.sym('interrmIn');
    interrTs  = SX.sym('interrTs');
    dummy_m_in_SP = SX.sym('dummy_m_in_SP');
    dummy_zV = SX.sym('dummy_zV');
    interrpS = SX.sym('interrpS');
    interrzV = SX.sym('interrzV');

    % Set states  x: differential, z: algebraic
    x = [];            
    for k=1:n
        x = [x; M(k); H(k)];
    end
    x = [x; Ms; Hs; MT; HT; w; Mb];
    % Control
    x = [x; interrmIn;interrTs;interrpS;interrzV];


    z = [m_in];
    for k=1:n
        z = [z; m(k);T(k);h(k);Tg(k);mg(k);Q(k);p(k);rho(k);beta(k);TB(k)];
    end

    z = [z; m_s; Ts; TBs; pS; hS; beta_s; rho_s; m_T; TT; TBT; pT; hT; ...
        beta_T; rho_T; W; P; m_U; TU; TU_real; TBU; hU; beta_U; Tc; m_c; ...
        hc; Qc];
 
    % Control
    z = [z; dummy_zm; dummy_m_in_SP;dummy_zV];
 
    % Initialize equatuions
    Alg  = [];
    diff = [];

    % Set function output variables
    x_keys = x;    
    z_keys = z;

    %% Additional inlet equations
    h_in  = cpW*(Tp-T_Ref);        % Note: Assuming input is liquid state

    % TC: PI Controller (Master)
    errTs = Ts_SP - Ts;                       % e: ySP-y
    m_in_SP = m_in_SP_0 + Kc_Ts*errTs + Kc_Ts/tauI_Ts*interrTs;

    % FC: I Controller
    errmIn = m_in_SP - m_in;                  % e: ySP-y
    z_m = z_m_0 + KI_m_in*interrmIn;          % u = u0 + KI*integral(e) 
    z_m_appl = fmax(0, fmin(1, z_m));         % u_bar = max(0, min(u, 1))

    init1 = m_in - z_m_appl*k_m*(p_pump-p(1));  
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

    U = U_Calc(beta(1), Tg(1), A_corr); % Pass numeric A_corr
    UA = A*U;
    
    VL    = (1-beta(1))*M(1)/rho(1);
    Pres1 = M(1)*R*T(1)/(MW*V) - p(1);
    Pres2 = p(1)*(V-VL) - beta(1)*M(1)*R*T(1)/(MW);
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

        U = U_Calc(beta(k), Tg(k), A_corr); % Pass numeric A_corr
        UA = A*U;

        VL    = (1-beta(k))*M(k)/rho(k);
        Pres1 = M(k)*R*T(k)/(MW*V) - p(k);
        Pres2 = p(k)*(V-VL) - beta(k)*M(k)*R*T(k)/(MW);
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

    U = U_Calc(beta(n), Tg(n), A_corr); % Pass numeric A_corr
    UA = A*U;

    VL    = (1-beta(n))*M(n)/rho(n);
    Pres1 = M(n)*R*T(n)/(MW*V) - p(n);
    Pres2 = p(n)*(V-VL) - beta(n)*M(n)*R*T(n)/(MW);
    Pres3 = rho(n) - M(n)/V;

    % Equations
    dMdt = m(n-1)-m(n);
    dHdt = m(n-1)*h(n-1) - m(n)*h(n) + Q(n);

    alg1  = mg(n)*cpG*(T0g-Tg(n)) - Q(n);
    alg2  = UA*((T0g-T(n))/2+(Tg(n)-T(n-1))/2) - Q(n);
    alg3  = m(n) - Cvd*(p(n)-pS);
    alg4  = if_else(cond1, Enth1, if_else(cond2, Enth3, Enth2));
    alg5  = M(n)*h(n) - H(n);
    alg6  = mg(n) - mG;
    alg7  = if_else(cond1, Pres1, if_else(cond2, Pres3, Pres2));
    alg8  = p(n) - 1/(k_p*rho_Ref)*(rho(n)-rho_Ref) - p_Ref;
    alg9  = cpW*(TB(n)-T_Ref) + beta(n)*dHvap - h(n);
    alg10 = 10^(A_anto-(B_anto/(TB(n)+C_anto))) - p(n);

    Alg  = [Alg;alg1;alg2;alg3;alg4;alg5;alg6;alg7;alg8;alg9;alg10];
    diff = [diff;dMdt;dHdt];

    %% Holdup and valve eq. after Hex (s)
    dHvap = dHvap0 + (cpW-cpS)*(TB_Ref-TBs); 
    % Switch logic
    cond1 = beta_s >= 1;
    cond2 = beta_s <= 0;

    Enth1 = cpW*(TBs-T_Ref) + dHvap + cpS*(Ts-TBs) - hS; 
    Enth2 = Ts - TBs;                                          
    Enth3 = cpW*(Ts-T_Ref) - hS;                             

    VLs   = (1-beta_s)*Ms/rho_s;
    Pres1 = Ms*R*Ts/(MW*Vs) - pS;
    Pres2 = pS*(Vs-VLs) - beta_s*Ms*R*Ts/(MW);
    Pres3 = rho_s - Ms/Vs;

    % Equations
    errpS = pS_SP - pS;                               % e: ySP-y
    zV = zV_0 + Kc_pS*errpS + Kc_pS/tauI_pS*interrpS + 1/tauT_pS*interrzV; 
    zV_appl = fmax(0, fmin(1, zV));                   % u_bar = max(0, min(u, 1))
    errzV = zV_appl - zV;

    dMsdt = m(n) - m_s;
    dHsdt = m(n)*h(n) - m_s*hS;

    alg1  = m_s - zV_appl*Kv*(pS-pT);
    alg2  = if_else(cond1, Enth1, if_else(cond2, Enth3, Enth2));
    alg3  = Ms*hS - Hs;
    alg4  = if_else(cond1, Pres1, if_else(cond2, Pres3, Pres2));
    alg5  = pS - 1/(k_p*rho_Ref)*(rho_s-rho_Ref) - p_Ref;
    alg6  = cpW*(TBs-T_Ref) + beta_s*dHvap - hS;
    alg7  = 10^(A_anto-(B_anto/(TBs+C_anto))) - pS;

    Alg  = [Alg;alg1;alg2;alg3;alg4;alg5;alg6;alg7];
    diff = [diff;dMsdt;dHsdt];

    %% Holdup before turbine (T)
    dHvap = dHvap0 + (cpW-cpS)*(TB_Ref-TBT); 
    % Switch logic
    cond1 = beta_s >= 1;
    cond2 = beta_s <= 0;

    Enth1 = cpW*(TBT-T_Ref) + dHvap + cpS*(TT-TBT) - hT; 
    Enth2 = TT - TBT;                                          
    Enth3 = cpW*(TT-T_Ref) - hT;                             

    VLT   = (1-beta_T)*MT/rho_T;
    Pres1 = MT*R*TT/(MW*VT) - pT;
    Pres2 = pT*(VT-VLT) - beta_T*MT*R*TT/(MW);
    Pres3 = rho_T - MT/VT;

    % Equations
    dMTdt = m_s - m_T;
    dHTdt = m_s*hS - m_T*hT;

    alg1  = if_else(cond1, Enth1, if_else(cond2, Enth3, Enth2));
    alg2  = MT*hT - HT;
    alg3  = if_else(cond1, Pres1, if_else(cond2, Pres3, Pres2));
    alg4  = pT - 1/(k_p*rho_Ref)*(rho_T-rho_Ref) - p_Ref;
    alg5  = cpW*(TBT-T_Ref) + beta_T*dHvap - hT;
    alg6  = 10^(A_anto-(B_anto/(TBT+C_anto))) - pT;

    Alg  = [Alg;alg1;alg2;alg3;alg4;alg5;alg6];
    diff = [diff;dMTdt;dHTdt];
    

    %% Turbine
    dwdt      = 1/Mg*(P-L-Dg*(w-w0));

    TurbAlg1  = m_T*sqrt(TT) - phiD*pT;
    TurbAlg2  = TU - TT*(pc/pT)^((R*10^5)/(cpS*MW*10^3));
    TurbAlg3  = W + eta*m_T*cpS*(TT-TU);
    TurbAlg4  = P + W;

    Alg  = [Alg; TurbAlg1; TurbAlg2; TurbAlg3; TurbAlg4];
    diff = [diff; dwdt];

    %% Condenser inlet flow (U)
    dHvap = dHvap0 + (cpW-cpS)*(TB_Ref - TBU);
    InletAlg1 = m_U - m_T;
    InletAlg2 = cpW*(TBU - T_Ref) + dHvap + cpS*(TU - TBU) - hU;
    InletAlg3 = cpW*(TBU - T_Ref) + beta_U*dHvap - hU;
    InletAlg4 = 10^(A_anto_c-(B_anto_c/(TBU+C_anto_c))) - pc;
    InletAlg5 = TU_real - TBU;

    Alg = [Alg; InletAlg1; InletAlg2; InletAlg3; InletAlg4; InletAlg5];


    %% Condenser (c)
    CondAlg1 = m_U - m_c;
    CondAlg2 = Qc - m_U*hU + m_c*hc;
    CondAlg3 = hc - cpW*(Tc-T_Ref);
    CondAlg4 = Tc - TBU;
    
    Alg = [Alg; CondAlg1;CondAlg2;CondAlg3;CondAlg4];


    %% Buffer tank
    dMbdt = m_c - m_in;

    diff = [diff; dMbdt];

    %% Append controller eq. 
    Dummy1 = dummy_zm - z_m_appl;
    Dummy2 = dummy_m_in_SP - m_in_SP;
    Dummy3 = dummy_zV - zV;
    Alg = [Alg; Dummy1; Dummy2; Dummy3];

    % Add error to ode
    diff = [diff; errmIn; errTs; errpS; errzV]; 

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