% Code to optimize value over Supercharge Pressure & tank temperature subcooling

global opts
dt = 0.1;

t0=0;                   %initial time of ignition
t_burn = 125;             %final time
tf=t_burn+t0;           %time when propelant is compeltely burned
t_range=[t0 tf];       %integration interval


N_opti = 10;
T_tank_range = 273.15 + linspace(-10,20,N_opti);
P_super_range = linspace(50e5,60e5,N_opti);

Heatmap = zeros(N_opti,N_opti);

for i=1:N_opti
    for j=1:N_opti
    
        T_init_ext = 283.15;    %K
        T_init_tank = T_tank_range(i);          %K
        opts.P_N2_init = P_super_range(i);      %Pa
        opts.V_N2_init = 0.05*opts.V_tank;      %m^3
        opts.T_N2_init = T_tank_range(i);       %K
        
        rho_liq = py.CoolProp.CoolProp.PropsSI('D','T',T_init_tank,'Q', 0,'NitrousOxide');
        rho_vap = py.CoolProp.CoolProp.PropsSI('D','T',T_init_tank,'Q', 1,'NitrousOxide');
        u_liq = py.CoolProp.CoolProp.PropsSI('U','T',T_init_tank,'Q', 0,'NitrousOxide');
        u_vap = py.CoolProp.CoolProp.PropsSI('U','T',T_init_tank,'Q', 1,'NitrousOxide');

        Fr_liq = opts.filling_ratio;        %WARNING : this is a volumic ratio and MUST NOT be confused with x, vapor massic ratio
        V_tank = opts.V_tank;

        m_liq = Fr_liq*V_tank*rho_liq;
        m_vap = (1-Fr_liq)*V_tank*rho_vap;

        m_ox_init=m_liq+m_vap;                                       %initial mass of oxidyzer in the tank
        U_total_init=m_liq*u_liq+m_vap*u_vap;                           %initial energy in the tank (both gaz and liquid phase)
        T_wall_init=T_init_ext;                                             %initial temperature of the tank wall (K)
        r_comb_chamber_init=opts.r_fuel_init;                            %initial radius in the combustion chamber (m)
        P_cc_init = 30*opts.P_atm_sl;                                   %initial pressure in combustion chamber (Pa)
        x_init=0;
        y_init=0;
        dxdt_init=0;
        dydt_init=0;
        r_throat_init = opts.D_throat/2;

        Initial_conditions=[m_ox_init; U_total_init; T_wall_init; r_comb_chamber_init; r_throat_init; P_cc_init; x_init; y_init; dxdt_init; dydt_init];%initial vector

        [t,state] = ode23s(@System_equations,t_range,Initial_conditions);%state1=m_tank_total, state2=U_tank_total,state3=T_tank_wall


        Heatmap(i,j) = value;
    end
end
