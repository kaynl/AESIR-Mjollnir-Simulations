tic

global opts

%% Set combustion parameters.
t0 = 0;                     % Initial time of ignition.
t_burn = 80;                % Final time.
tf = t_burn + t0;           % Time when propelant is completely burned.
t_range = [t0 tf];          % Integration interval.


% TODO: Figure out a way to run this script with a given choice of parametr
%       values instead of having to change it in the file every time.
% TODO: Remove the pressurization system.



%% Initialization.

disp("---------------------------------")
disp("Intitialization...") 
disp("---------------------------------")
disp(" ")

% TODO: move settings to main interface.
[T_ext_model, ~, ~, ~] = atmoscoesa(0);
T_init_ext = 283;                                   % Initial external temperature (K).
opts.dT_ext = T_ext_model - T_init_ext;             % Difference between the initial temperature and the model temperature (K).
T_init_tank = 286;                                  % Initial tank temperature is room temperature (K).

rho_liq = py.CoolProp.CoolProp.PropsSI('D','T',T_init_tank,'Q', 0,'NitrousOxide');      % Density for liquid N2O (kg/m^3).
rho_vap = py.CoolProp.CoolProp.PropsSI('D','T',T_init_tank,'Q', 1,'NitrousOxide');      % Density for vapor N2O (kg/m^3).
u_liq = py.CoolProp.CoolProp.PropsSI('U','T',T_init_tank,'Q', 0,'NitrousOxide');        % Specific internal energy for liquid N2O (J/kg).
u_vap = py.CoolProp.CoolProp.PropsSI('U','T',T_init_tank,'Q', 1,'NitrousOxide');        % Specific internal energy for vapor N2O (J/kg).

fill_ratio = opts.filling_ratio;        % WARNING : this is a volumic ratio and MUST NOT be confused with x, massic vapor ratio.
V_tank = opts.V_tank;                               % The volume of the tank (m^3).

m_liq = fill_ratio * V_tank * rho_liq;              % The liquid mass is the liquid volume in the tank times the liquid density (kg).
m_vap = (1 - fill_ratio) * V_tank * rho_vap;        % The liquid mass is the remaining volume in the tank times the vapor density (kg).

m_ox_init = m_liq + m_vap;                          % The initial mass of the oxidizer in the tank is the sum of liquid and vapor mass (kg).
U_total_init = m_liq * u_liq + m_vap * u_vap;       % The initial energy in the tank is the sum of liquid and vapor mass times energy (J).
T_wall_init = T_init_ext;                           % The initial temperature of the tank wall is equal to the external temperature (K).
r_cc_init = opts.r_fuel_init;                       % The initial radius of the combustion chamber is equal to the initial radius of the fuel port (m).
% TODO: why 20?
P_cc_init = 15 * opts.P_atm_sl;                     % The initial pressure in the combustion chamber is twenty times the atmospheric pressure (Pa).
% TODO: what is D_throat, and why is the initial radius of the nozzle throat half of this quantity?
r_throat_init = opts.D_throat / 2;                  % The initial radius of the nozzle throat is...

x_init = 0;
y_init = 0;
dxdt_init = 0;
dydt_init = 0;

%% Solve differential equations.

disp("---------------------------------")
disp("Solving differential equations...") 
disp("---------------------------------")
disp(" ")

% Initialization.
% tol = odeset('RelTol',1e-5,'AbsTol',1e-6);
initial_conditions = [m_ox_init; U_total_init; T_wall_init; r_cc_init; r_throat_init; P_cc_init; x_init; y_init; dxdt_init; dydt_init];

% Solve initial value problem for ODE.
opts.tank_state = 100;
% ode_opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-8);
% ode_opts = odeset('MaxStep', 0.05);
% [t, state] = ode23tb(@System_equations, t_range, initial_conditions, ode_opts);
[t, state] = ode23tb(@System_equations, t_range, initial_conditions);

% Retrieve results.
m_ox_total = state(:, 1);                                                           % The total mass in the tank at each time step.
empty_tank = find(m_ox_total < 0);
empty_tank_mask = ones(length(t), 1);
empty_tank_mask(empty_tank) = 0;
m_ox_total = m_ox_total .* empty_tank_mask;                                         % Set negative mass to zero.

U_tank_total = state(:, 2) .* empty_tank_mask;                                      % The total energy inside the tank at each time step.
T_tank_wall = state(:, 3);                                                          % The wall temperature in the tank at each time step.
r_cc = state(:, 4);                                                                 % The radius of the combustion chamber at each time step.
r_throat = state(:, 5);
P_cc = state(:, 6) .* empty_tank_mask + opts.P_atm_sl .* (1 - empty_tank_mask);     % The pressure in the combustion chamber at each time step.
x = state(:, 7);                                                                    % The position on the x-axis (m).
y = state(:, 8);                                                                    % The position on the y-axis (m).
vx = state(:, 9);                                                                   % The velocity along the x-axis (m/s).
vy = state(:, 10);                                                                  % The velocity along the y-axis (m/s).

save('simulation_pre_compute')

%% Post-compute recuperation of data.

disp("---------------------------------")
disp("Post-compute calculations...") 
disp("---------------------------------")
disp(" ")

% Post recuperation of data.
N = length(vx);                     % Save every tenth data point.
T_tank = zeros(1, N);               % The temperature inside the tank (K).
V_tank = zeros(1, N);               % The volume of the tank (m^3).
V_liq = zeros(1, N);                % The volume of the liquid in the tank (m^3).
P_tank = zeros(1, N);               % The pressure inside the tank (Pa).
x_vap = zeros(1, N);                % The massic vapor ratio in the tank.
mf_ox = zeros(1, N);                % The oxidizer mass flow (kg/s).
mf_fuel = zeros(1, N);              % The fuel mass flow (kg/s).
OF = zeros(1, N);                   % The oxidizer-fuel ratio.
mf_throat = zeros(1, N);            % The mass flow through the throat.
T_cc = zeros(1, N);                 % The temperature inside the combustion chamber (K).
Tr = zeros(1, N);                   % The thrust (N).
ax = zeros(1, N);                   % The acceleration along the x-axis (m/s^2).
ay = zeros(1, N);                   % The acceleration along the y-axis (m/s^2).
Me = zeros(1, N);                   % The exhaust Mach.
Ve = zeros(1, N);                   % The exhaust velocity (m/s).
Pe = zeros(1, N);                   % The exhaust pressure (Pa).
At = pi * r_throat.^2;              % The throat area (m^2).
cp_air = zeros(1, N);               % The specific heat of air (J/kg/K).
P_N2O = zeros(1, N);                % The saturation pressure of N20 (Pa).
m_fuel = zeros(1, N);               % The fuel mass (kg).
h_liq = zeros(1, N);                % The heat transfer coefficient for the liquid N2O (W/m^2/K).
h_gas = zeros(1, N);                % The heat transfer coefficient for the gaseous N2O (W/m^2/K).
h_air_ext = zeros(1, N);            % The heat transfer coefficient for air (W/m^2/K).
m_tot = zeros(1, N);                % The total mass (kg).

% TODO: I will figure out the rest of the code later.

% Get the temperature (K), speed of sound (m/s), pressure (Pa), and density (kg/m^3) at each height y_i.
[T_ext, speed_of_sound, P_ext, rho_ext] = atmoscoesa(y);

A_port = pi * r_cc.^2;      % The area of the combustion chamber port (m^2).

gam = opts.gamma_combustion_products;               % The heat capacity ratio.
Mw = opts.molecular_weight_combustion_products;     % The molecular weight of products (kg/mol).
R = opts.R;                                         % The universal gas constant (J/K/mol).        

% pause(10) % Why would we?

opts.tank_state = 100;
for i = 1:length(T_tank)
    [dstatedt, m_tot(i), h_liq(i), h_gas(i), h_air_ext(i)] = System_equations(t(i), state(i, :));
    ax(i) = dstatedt(9);
    ay(i) = dstatedt(10);
    
    if ~ismember(i, empty_tank)
        T_tank(i) = Tank_Temperature(U_tank_total(i), m_ox_total(i));
        x_vap(i) = x_vapor(U_tank_total(i), m_ox_total(i), T_tank(i));
        
        rho_vap = fnval(opts.RhoG_T_N2O_spline, T_tank(i));
        rho_liq = fnval(opts.RhoL_T_N2O_spline, T_tank(i));
        % TODO: rename to remaining_liquid?
        Tank_state = (1 - x_vap(i)) * m_ox_total(i) / (rho_liq * opts.V_tank) * 100;
        % TODO: pressurization.
        T_N2 = T_tank(i);
        % TODO: pressurization.
        V_N2 = x_vap(i) * m_ox_total(i) / rho_vap;
        V_tank(i) = m_ox_total(i) * (x_vap(i) / rho_vap + (1 - x_vap(i)) / rho_liq);    % Compute the total remaining volume in the tank, recall that volume = mass / density.
        V_liq(i) = m_ox_total(i) * (1 - x_vap(i)) / rho_liq;                            % Compute the remaining liquid volume in the tank.
        
        % TODO: why is this 10^6 here, and 10^5 in other places?
        P_N2O(i) = fnval(opts.Psat_N2O_spline, T_tank(i)) * 10^6;                                                               % Compute the 
        P_tank(i) = P_N2O(i);
        
        mf_ox(i) = Mass_flow_oxidizer(T_tank(i), P_tank(i), P_cc(i));
        G_Ox = mf_ox(i)/ A_port(i);
        mf_fuel(i) = Mass_flow_fuel(G_Ox, r_cc(i));
        OF(i) = mf_ox(i)/ mf_fuel(i);
        mf_throat(i) = Mass_flow_throat(P_cc(i), OF(i), At(i)); 
        c_star = interp1q(opts.OF_set, opts.C_star_set, OF(i));
        T_cc(i) = gam * (2 / (gam + 1))^((gam + 1) / (gam - 1)) * c_star^2 * Mw / R;
        Me(i) = ExhaustMach(opts, At(i));
        Pe(i) = ExhaustPressure(P_cc(i), P_ext(i), Me(i), opts);
        Ve(i) = ExhaustSpeed(T_cc(i), Pe(i), P_cc(i), opts);
        Tr(i) = opts.combustion_efficiency * Thrust(mf_throat(i), Ve(i), P_ext(i), Pe(i), opts);
        cp_air(i) = py.CoolProp.CoolProp.PropsSI('C', 'P', P_ext(i), 'T', (T_ext(i) + T_tank_wall(i)) / 2, 'Air');
        m_fuel(i) = opts.rho_fuel * opts.L_fuel * pi * (opts.D_cc_int^2 / 4 - r_cc(i)^2);
    else
        % Tank is empty.
        index = empty_tank(1) - 1;          % This is the last observation with non-empty tank.
        T_tank(i) = T_tank(index);          % Take the temperature inside the tank to be equal to the final temperature.
        P_tank(i) = opts.P_atm_sl;          % Take the pressure inside the tank to be equal to the atmospheric pressure.
        T_cc(i) = T_cc(index);              % Take the temperature inside the combustion chamber to be equal to the final temperature.
    end
end

save('simulation_results')

toc