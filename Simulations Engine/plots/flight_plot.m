%% Plots of combustion part

global opts

lignes=2;
colonnes=4;
figure(1)
sgtitle("Flight Parameters",'FontSize', 20, 'Color','Red','FontWeight','bold')

m_fuel = opts.rho_fuel*opts.L_fuel*pi*(opts.D_cc_int.^2/4-r_cc.^2);


v=sqrt(vx.^2+vy.^2);
acc = sqrt(ax.^2+ay.^2);

Cd = drag_coefficient_model(v, speed_of_sound);
D= Cd .* 0.5 .* rho_ext .* v.^2 .* opts.surface;

Dx = D.*vx./v;
Dy = D.*vy./v;

%% Speed
subplot(lignes,colonnes,1)
plot(t,sqrt(vx.^2+vy.^2),t,vx,t,vy)
title("Speed Over Time")
xlabel("Time (s)")
ylabel("Speed (m/s)")
legend("Speed","V_x","V_y")

%% Drag
subplot(lignes,colonnes,2)
plot(t,D,t,Dx,t,Dy)
title("Drag Over Time")
xlabel("Time (s)")
ylabel("Drag (N)")
legend("Total Drag","D_x","D_y")

%% Acceleration
subplot(lignes,colonnes,3);
plot(t,acc/opts.g,t,ax/opts.g,t,ay/opts.g)
title("Acceleration Over Time")
xlabel("Time (s)")
ylabel("Acceleration (g)")
legend("Total Acceleration","a_x","a_y")

%% Mach Number

subplot(lignes,colonnes,4);
plot(t,v./speed_of_sound)
title("Mach Over Time")
xlabel("Time (s)")
ylabel("Mach")

%% Trajectory
subplot(lignes,colonnes,5)
plot(x/1000,y/1000,x/1000,y*0/1000+opts.design_altitude/1000,'--',x/1000,y*0/1000+opts.required_altitude/1000)
title("Trajectory")
xlabel("Longitudinal distance (km)")
ylabel("Height (km)")
legend("Trajectory","14 km","12 km")

%% Height
subplot(lignes,colonnes,6)
plot(t,y/1000,t,y*0/1000+opts.design_altitude/1000,'--',t,y*0/1000+opts.required_altitude/1000)
title("Heigth Over Time")
xlabel("Time (s)")
ylabel("Height (km)")
legend("Trajectory","14 km","12 km")

%% Air Pressure Over Time

subplot(lignes,colonnes,7);
plot(t, P_ext/10^5)
title("Air Pressure Over Time")
xlabel("Time (s)")
ylabel("Air Pressure (bar)")

%% Dynamic Pressure

subplot(lignes,colonnes,8);
plot(y/1000, 0.5*rho_ext.*v.^2./10^5)
title("Dynamic Pressure Over Altitude")
xlabel("Altitude (km)")
ylabel("Dynamic Pressure (bar)")

