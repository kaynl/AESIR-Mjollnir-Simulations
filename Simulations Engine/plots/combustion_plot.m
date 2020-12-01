%% Plots of combustion part

global opts

lignes=2;
colonnes=4;
figure(2)
sgtitle("Combustion Parameters",'FontSize', 20, 'Color','blue','FontWeight','bold')

t_burn = 15;
t_comb=t(find(t<t_burn));

%% Radius Port Over time
subplot(lignes,colonnes,1);
plot(t_comb,r_cc(find(t<t_burn))*1000,t_comb,0*r_cc(find(t<t_burn))+opts.D_cc_int/2*1000,t_comb,0*r_cc(find(t<t_burn))+opts.fuel_margin_radius*1000)
title("Port Radius Over Time")
xlabel("Time (s)")
ylabel("Port Radius (mm)")
lgd = legend("Port Radius","External Tank Limit","Fuel Margin Limit");
lgd.Location = 'southeast';

%% Tank Temperature Over Time
subplot(lignes,colonnes,2)

v=sqrt(vx.^2+vy.^2)
T_wall_ext = T_ext + v.^2./(2.*cp_air);

plot(t_comb,T_tank(find(t<t_burn)), t_comb, T_tank_wall(find(t<t_burn)), t_comb, T_ext(find(t<t_burn)))
title("Tank Temperatures Over Time")
xlabel("Time (s)")
ylabel("Temperature Tank (K)")

lgd = legend("Tank Int Temperature","Tank Ext Wall Temperature","Ext Air Temperature");
lgd.Location = 'southwest';

%% Frictional Temperature Over Time
subplot(lignes,colonnes,3)
plot(t_comb, T_wall_ext(find(t<t_burn)))
title("Frictional Temperature Over Time")
xlabel("Time (s)")
ylabel("Friction Temperature (K)")

%% Fuel Regression Rate
subplot(lignes,colonnes,4)
a = opts.reg_a;
n = opts.reg_n;
A_fuel = pi*r_cc'.^2;
G_Ox = mf_ox./A_fuel;

drdt=a*(G_Ox).^n;
plot(t_comb,drdt(find(t<t_burn))*1000)
title("Fuel Regression Rate Over Time")
xlabel("Time (s)")
ylabel("Regression Rate (mm/s)")

%% Mass flow Ox + Fuel
subplot(lignes,colonnes,5)
plot(t_comb,mf_ox(find(t<t_burn)),t_comb,mf_fuel(find(t<t_burn)),t_comb,mf_throat(find(t<t_burn)))
title("Mass Flows Over Time")
xlabel("Time (s)")
ylabel("Mass Flow (kg/s)")
lgd = legend("Oxidizer","Fuel","Throat");
lgd.Location = 'east';

%% O/F Ratio

subplot(lignes,colonnes,6)
plot(t_comb,OF(find(t<t_burn)))
title("O/F Ratio Over Time")
xlabel("Time (s)")
ylabel("O/F")

%% Pressure in CC and Tank over Time


subplot(lignes,colonnes,7)
plot(t_comb,P_tank(find(t<t_burn))/10^6,t_comb,P_cc(find(t<t_burn))/10^6,t_comb,Pe(find(t<t_burn))/10^6)
title("Pressure Over Time")
xlabel("Time (s)")
ylabel("Pressure (MPa)")
lgd = legend("Tank Pressure","CC Pressure", "Exhaust Pressure");
lgd.Location = 'southwest';

%% Tank fillness

subplot(lignes,colonnes,8)
plot(t_comb,m_ox_total(find(t<t_burn))/opts.rho_ox*1000,t_comb,opts.V_tank*1000+0*m_ox_total(find(t<t_burn))/opts.rho_ox,0,0)
title("NO2 Volume Over Time")
xlabel("Time (s)")
ylabel("Tank Volume (L)")
lgd = legend("Tank Volume","Full Tank Volume");
lgd.Location = 'east';

