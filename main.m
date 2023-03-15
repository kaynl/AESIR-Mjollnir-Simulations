%% This file is the main user interface.
setup;

%% User settings.
simulate = false;             % True if the simulation should be run, if false it will load the most recent simulation.
process_data = false;        % TODO: it would be nice to integrate this more properly into the main.
data_name = "Datasets/HT-2/ht21_large.mat";

plot_data = true;           % True if the data should be plot together with the simulations.

%% Model parameters.
full_duration = false;       % True if the tank parameters should be set to a full-duration burn, otherwise short-duration parameters are used.
Cd = 0.66;                  % Discharge coefficient.
dr_thdt = 0.1e-3;           % Constant approximation of regression rate (m/s).

%% Vehicle parameters.
n_inj = 80;                 % Number of injector holes.

%% Set options.
set_options(plot_data, full_duration, Cd, dr_thdt, n_inj)

%% Run or load simulations.
if simulate
    % Run simulations.
    simulate_engine;
end
simulation = load("simulation_results.mat");

%% Data processing.
if process_data
    disp("")          %% TODO: data processing.
else
    data = load(data_name); 
    data = data.data;
end

%% Plot the simulation results.
plot_results;
