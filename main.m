%% This file is the main user interface.

%% User settings.
simulate = true;            % True if the simulation should be run, if false it will load the most recent simulation.
process_data = true;        % TODO: it would be nice to integrate this more properly into the main.

plot_data = true;           % True if the data should be plot together with the simulations.

%% Model parameters.
full_duration = true;       % True if the tank parameters should be set to a full-duration burn, otherwise short-duration parameters are used.
Cd = 0.66;                  % Discharge coefficient.
dr_thdt = 0.1e-3;           % Constant approximation of regression rate (m/s).

%% Vehicle parameters.
n_inj = 80;                 % Number of injector holes.

%% Run setup and set manual parameters.
setup;
global opts
opts.plot_data = plot_data;
opts.full_duration = full_duration;
opts.Cd = Cd;
opts.dr_thdt = dr_thdt;
opts.n_inj = n_inj;

if full_duration
    opts.L_tank = 1.83;     % Tank length for full-duration burn (m).
else
    opts.L_tank = 0.73;     % Tank length for short-duration burn (m).
end

%% Run or load simulations.
if simulate
    % Run simulations.
    Main_engine;
else
    % Load latest result.
    load("simulation_results.mat");
end

%% Data processing.
if process_data
    disp()          %% TODO: data processing.
end

%% Plot the simulation results.
disp("---------------------------------")
disp("Plotting...") 
disp("---------------------------------")
disp(" ")

plot(plot_data);
