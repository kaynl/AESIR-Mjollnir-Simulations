%% Thrust plots.

global opts

rows = 2;
columns = 3;
figure(2)
sgtitle("Thrust", 'FontSize', 20, 'Color', 'Green', 'FontWeight', 'bold')

if opts.full_duration
    t_burn = 25;
else
    t_burn = 10;
end
sim_ind = find(simulation.t < t_burn);
t_sim = simulation.t(sim_ind);
t_data = linspace(0, t_burn, 1000);

%% Thrust.
subplot(rows, columns, 1);
plot(t_sim, simulation.Tr(sim_ind) / 1000)
hold on
plot(t_sim, opts.combustion_efficiency * simulation.Tr(sim_ind) / 1000)
if opts.plot_data
    plot(t_data, data.THRUST_I(t_data), '--', 'Color', '#EDB120');
end
hold off

title("Thrust")
xlabel("Time (s)")
ylabel("Thrust (kN)")
if opts.plot_data
    legend("Ideal", "Estimated", "Actual", 'Location', 'southeast')
else
    legend("Ideal", "Estimated", 'Location', 'southeast')
end
axis padded

%% Specific impulse.
subplot(rows, columns, 2)
Isp = simulation.Tr ./ (opts.g .* simulation.mf_throat);
plot(t_sim, Isp(sim_ind))

title("Specific impulse")
xlabel("Time (s)")
ylabel("Isp (s)")
axis padded

%% CC temperature.
subplot(rows, columns, 3)
plot(t_sim, simulation.T_cc(sim_ind))

title("CC temperature")
xlabel("Time (s)")
ylabel("Temperature (K)")
axis padded

%% Exhaust speed.
subplot(rows, columns, 4)
plot(t_sim, simulation.Ve(sim_ind))

title("Exhaust speed")
xlabel("Time (s)")
ylabel("Speed (m/s)")
axis padded

%% Nozzle exit area.
subplot(rows, columns, 5)
plot(t_sim, simulation.At(sim_ind))

title("Nozlle exit area")
xlabel("Time (s)")
ylabel("Area (m²)")
axis padded

%% Exhaust mach.
subplot(rows, columns, 6)
plot(t_sim, simulation.Me(sim_ind))

title("Exhaust mach")
xlabel("Time (s)")
ylabel("Mach")
axis padded
