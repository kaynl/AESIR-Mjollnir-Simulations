%% Plots of thrust part

global opts

lignes=2;
colonnes=3;
figure(3)
sgtitle("Thrust Parameters",'FontSize', 20, 'Color','Green','FontWeight','bold')

t_burn = 15;
t_comb=t(find(t<t_burn));

%% Thrust Over Time

subplot(lignes,colonnes,1);
plot(t_comb,Tr(find(t<t_burn))/1000, t_comb, opts.combustion_efficiency*Tr(find(t<t_burn))/1000)
title("Thrust Over Time")
xlabel("Time (s)")
ylabel("Thrust (kN)")
lgd = legend("Ideal Thrust","Real Thrust Estimation")
lgd.Location = 'southwest';

%% Isp Over Time

Isp = Tr./(opts.g.*mf_throat);
subplot(lignes,colonnes,2)
plot(t_comb,Isp(find(t<t_burn)))
title("Isp Over Time")
xlabel("Time (s)")
ylabel("Isp (s)")




%% CC Temperature Over Time
subplot(lignes,colonnes,3)
plot(t_comb,T_cc(find(t<t_burn)))
title("CC Temperature Over Time")
xlabel("Time (s)")
ylabel("CC Temperature (K)")


%% Exhaust Speed Over Time

subplot(lignes,colonnes,4)
plot(t_comb,Ve(find(t<t_burn)))
title("Exhaust Speed Over Time")
xlabel("Time (s)")
ylabel("Exhaust Speed (m/s)")

%% Nozzle Exit Area Over Time

subplot(lignes,colonnes,5)
plot(t_comb,At(find(t<t_burn)))
title("Nozlle Exit Area Over Time")
xlabel("Time (s)")
ylabel("Nozzle Exit Area (m²)")

%% Exhaust Mach Over Time

subplot(lignes,colonnes,6)
plot(t_comb,Me(find(t<t_burn)))
title("Exhaust Mach Over Time")
xlabel("Time (s)")
ylabel("Exhaust Mach")
