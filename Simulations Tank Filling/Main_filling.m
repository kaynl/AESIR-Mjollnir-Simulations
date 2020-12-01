run('./../setup')

global opts
global dt

dt = 0.001;

t0=0;                           %initial time of filling
t_filling =10*60;               %final time
tf=t_filling+t0;                %arbitrary time 
t_range=t0:0.01:tf;                %integration interval

%Initial conditions
disp("-----------------------")
disp("Intitialization") 
disp("-----------------------")
disp(" ")

T_init = opts.T_ext;    %K
V_tank = opts.V_tank;   %m^3

rho_vap = py.CoolProp.CoolProp.PropsSI('D','P',opts.P_storage_tank, 'Q', 1,'NitrousOxide');

% V_liq_init = 0;
V_vap_init = V_tank;
m_liq_init = 0;
m_vap_init = opts.P_storage_tank*V_tank/(opts.r_ox*T_init);%rho_vap*V_vap_init;%

u=py.CoolProp.CoolProp.PropsSI('U','T',T_init,'Q', 1,'NitrousOxide');
U_tot_init = m_vap_init*u;

Initial_conditions=[m_liq_init; m_vap_init; U_tot_init];   %initial vector

%Solve initial value problem for ODE
disp("-----------------------")
disp("Solving Differential Eq") 
disp("-----------------------")
disp(" ")
[t,state] = ode45(@System_of_equations_filling,t_range,Initial_conditions);


m_liq = state(:,1);
m_vap = state(:,2);
U_tot = state(:,3);

%%

Results_plots;