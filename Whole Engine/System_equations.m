function [state_vector, m_tot, h_liq, h_gas, h_air_ext]  = System_equations(t,u)

global opts


% TODO: Figure out some way of selecting which parts of the model to use.


m_ox = u(1,:);
U_total = u(2,:);
T_wall = u(3,:);
r_cc = u(4,:);
r_throat = u(5,:);
P_cc = u(6,:);
x = u(7,:);
y = u(8,:);
dxdt = u(9,:);
dydt = u(10,:);


At = pi*r_throat.^2;
dr_thdt = 0.1e-3;   %0.1mm/s regression rate : constant approximation, TO BE CORRECTED AFTER STATIC FIRE
V_rocket = sqrt(dxdt^2+dydt^2);

[T_ext_model, speed_of_sound, P_ext, rho_ext] = atmoscoesa(y);
T_ext = T_ext_model - opts.dT_ext;
flight_state = opts.flight_state;       %0 for static test and 1 for in flight test


a = opts.reg_a;
n = opts.reg_n;

D_tank_ext = opts.D_ext_tank;
D_tank_int = opts.D_int_tank;
L_tank = opts.L_tank;
m_wall = opts.rho_alu*pi*(D_tank_ext^2-D_tank_int^2)/4*L_tank;
c_wall = opts.alu_thermal_capacity; 
L_cc = opts.L_cc;
L_fuel = opts.L_fuel;
rho_fuel = opts.rho_fuel;
D_cc_int = opts.D_cc_int;
gam = opts.gamma_combustion_products;
Mw = opts.Molecular_weigth_combustion_products;
R = opts.R;

%Tank Temperature Calculation
T_tank=Tank_Temperature(U_total,m_ox);

x=x_vapor(U_total,m_ox,T_tank); %x_vapor computed thanks to the internal tank temperature
% disp("x : "+x)
rho_liq = fnval(opts.RhoL_T_NO2_spline,T_tank);
rho_vap = fnval(opts.RhoG_T_NO2_spline,T_tank);

Tank_state = (1-x)*m_ox/(rho_liq*opts.V_tank)*100;

%Tank Pressure Calculation (Psat at T_int) with supercharge
T_N2 = T_tank;
V_vap = (x*m_ox/rho_vap);
V_N2 = V_vap;

% disp("V_N2 : "+V_N2*1000+" L")

P_N2 = (opts.P_N2_init*opts.V_N2_init/opts.T_N2_init)*T_N2/V_N2;

P_N2O = fnval(opts.Psat_NO2_spline,T_tank)*10^6;
% disp("P_N2O : "+P_N2O/10^5+" bars")


P_tank = P_N2O+abs(opts.P_N2_init-fnval(opts.Psat_NO2_spline,opts.T_N2_init)*10^6)*Tank_state/100;
% disp("P_tank : "+P_tank/10^5+" bars")

% disp("Pcc : "+P_cc/10^5+" bars")
 
%Massflows oxidizer/fuel/throat
mf_ox=Mass_flow_oxidizer(T_tank,P_tank,P_cc);%outlet mass flow

% disp("mf ox : "+mf_ox)
A_fuel = pi*r_cc^2;
G_Ox = mf_ox/A_fuel;
mf_fuel = Mass_flow_fuel(G_Ox,r_cc);
OF = mf_ox/mf_fuel;
% disp("OF : "+OF)

mf_throat = Mass_flow_throat(P_cc,OF,At);

%Heat Flux Interior of the tank
[Qdot_w_t, h_liq, h_gas]=HeatFlux_wall_tank(P_N2O, x, T_wall, T_tank);%Thermal heat flux from the wall to the tank

%Heat Flux Exterior of the tank
[Qdot_ext_w, h_air_ext] = HeatFlux_ext_wall(V_rocket,T_ext,P_ext,T_wall, flight_state,opts);

%Enthalpie outlet tank
h_outlet = py.CoolProp.CoolProp.PropsSI('H','T',T_tank,'Q', 0,'NitrousOxide');

%Combustion Chamber Volume calculation

V_cc = (L_cc*pi*D_cc_int^2/4)-(L_fuel*pi*(D_cc_int^2/4-r_cc^2));
m_fuel = rho_fuel*L_fuel*pi*(D_cc_int^2/4-r_cc^2);

c_star = interp1q(opts.OF_set,opts.C_star_set,OF);
RTcc_Mw = gam*(2/(gam+1))^((gam+1)/(gam-1))*c_star^2;
T_cc = RTcc_Mw*Mw/R;

%Equations that model the tank :
dmtotaldt=-mf_ox;
dUtotaldt=-mf_ox*h_outlet+Qdot_w_t;
% disp("Qdot_ext_w : "+Qdot_ext_w)
% disp("Qdot_w_t : "+Qdot_w_t)
dTwalldt=(Qdot_ext_w-Qdot_w_t)/(m_wall*c_wall);
drdt=a*(G_Ox)^n;
dP_ccdt=(mf_fuel+mf_ox-mf_throat)*RTcc_Mw/(V_cc);


%Exhaust parameters

Me = ExhaustMach(opts,At);
Pe = ExhaustPressure(P_cc,P_ext,Me,opts);
Ve = ExhaustSpeed(T_cc,Pe,P_cc,opts);
% disp("Exhaust Speed (m/s) : "+Ve)



if Tank_state<1
    % Eq of Motion
    m_tot = m_fuel+opts.dry_mass;
    [d2xdt2, d2ydt2] = EqofMotion(0,0,m_fuel,y,dxdt,dydt,speed_of_sound,rho_ext, flight_state,opts);
    if y<0
        dxdt = 0;
        dydt = 0;
        d2xdt2 = 0;
        d2ydt2 = 0;
    end
    state_vector = [0;0;0;0;0;0;dxdt;dydt;d2xdt2;d2ydt2];
else
    
    % Eq of Motion
    
%     Cs = ThrustCoefficient(Pe,P_cc,P_ext,At,opts);
%     Tr = opts.combustion_efficiency*Thrust_using_Cs(mf_throat,Cs,OF,opts);
    Tr = opts.combustion_efficiency*Thrust(mf_throat,Ve,P_ext,Pe,opts);
%     disp("Thrust (kN) : "+Tr/1000)
    m_tot = m_ox+m_fuel+opts.dry_mass;
    [d2xdt2, d2ydt2] = EqofMotion(Tr,m_ox,m_fuel,y,dxdt,dydt,speed_of_sound,rho_ext, flight_state,opts);
    
    
    state_vector=[dmtotaldt; dUtotaldt; dTwalldt; drdt; dr_thdt; dP_ccdt;dxdt;dydt;d2xdt2;d2ydt2];
end



% % Disp Section
% disp("dUdt : "+dUtotaldt)
% disp("dmdt :"+dmtotaldt)
% disp("dP (bars) : "+(P_tank-P_cc)/10^5)
% disp("dTwalldt : "+dTwalldt)
% disp("dPccdt (bar/s): "+dP_ccdt/10^5)
% disp("mass flow ox : "+mf_ox)
% disp("mf throat : "+mf_throat)
% disp("mf fuel : "+mf_fuel)
% disp("Thrust : "+Tr/1000+" kN")
% disp("Speed : "+sqrt(dxdt^2+dydt^2))
% disp("Acceleration : "+sqrt(d2xdt2^2+d2ydt2^2))
% disp("a_x (g) : "+d2xdt2./opts.g)
% disp("a_y (g) : "+d2ydt2./opts.g)
% disp("Height : "+y)
if Tank_state<opts.tank_state-5 %|| Tank_state<1.2
    opts.tank_state = opts.tank_state-5;
    disp("Tank State : "+round(Tank_state)+" % full")
    
end
% disp(" ")

end

