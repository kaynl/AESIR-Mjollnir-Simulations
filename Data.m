%%% This file gathers all engine data needed for the simulation

% TODO: Split data into parameters and "constants" based on how we need
%       to use them.
% TODO: Create a good way to run a specific simulation with a certain set
%       of parameter choices.

global opts

opts.flight_state = 1;          % 0 if test fire on the ground and 1 if launch.
opts.filling_ratio = 0.95;      % Tank filling ratio.
opts.launch_angle = 87;         % Self explanatory (in °).

opts.drag_coefficient = 0.5;     
opts.combustion_efficiency = 0.9;
opts.T_ext = 293;               % (K) Exterior temperature (20°C).

%% Physical constants.

opts.g = 9.81;                                  % Gravitational constant (m/s^2).
opts.R = 8.314;                                 % Universal gas constant (J/K/mol).

opts.P_atm_sl = 101325;                         % Atmospheric pressure (N/m2).
opts.stephan_cst = 5.67e-8;                     % Stephan-Boltzman constant (W/m2/K4).
opts.eber_parameter = 0.89;                     % Eber parameter for vertex angle between 20-50 degrees.
opts.Molecular_weight_air = 28.9647e-3;         % Molecular weight of air (kg/mol).
opts.r_air = opts.R / opts.Molecular_weight_air;


%% Requirements.

opts.design_altitude = 14000;          % Designed altitude to reach (m).
opts.required_altitude = 12000;        % Mission requirements (m).


%% Mass.

opts.parachute_mass = 10;
opts.electronics_mass = 2.3;
opts.bodyTube_mass = 7;
opts.payload_mass = 2;

opts.propulsionSystem = 24.504;

opts.dry_mass = opts.parachute_mass + opts.electronics_mass + opts.bodyTube_mass + opts.payload_mass + opts.propulsionSystem;

opts.m_ox_init = 24.5;          % Oxidizer mass (kg).
opts.m_fuel_init = 3.1;         % Fuel mass (kg).
opts.rho_ox = 785;              % Oxidizer density (kg/m^3).

%% Tank geometry

% TODO: Turn on/off big tank.

opts.D_ext_tank = 10e-2;    % Big: 16-2, small: 10e-2.                                     % Tank external diameter (m).
opts.e_tank = 3.5e-3;                                                   % Tank thickness.
opts.D_int_tank = opts.D_ext_tank - 2 * opts.e_tank;    %9.42e-2;       % Tank internal diameter (m).
opts.L_tank = 0.73;     % Big: 1.83, small: 0.73.                                      % Tank length (m).
opts.V_tank = pi * (opts.D_int_tank)^2 / 4 * opts.L_tank;   %33.1e-3;   % Tank volume (m^3) (present in Tank_Temperature_finder_fct).
opts.surface = pi * (opts.D_ext_tank)^2 / 4;                            % Rocket surface.

%% Kastrullen.

opts.L_kastrullen = 35e-2;  % Length of Kastrullen.

%% Injector geometry.
              
opts.r_inj = 1.2e-3 / 2;        % Injector radius (m).
opts.Cd = 0.66;                 % Discharge coefficient.
opts.L_inj = 15e-3;             % Injector plate thickness (m).
opts.n_inj = 80;                % Number of injectors.

opts.r_inj_plate = 30e-3;       % m
opts.mass_inj = 0.271;          % kg
opts.e_inj = 0.013;             % m

%% Combustion chamber geometry.

opts.D_cc_ext = 15.2e-2;                        % Combustion chamber external diameter (m).
opts.e_cc = 4e-3;
opts.D_cc_int = opts.D_cc_ext-2 * opts.e_cc;    % Combustion chamber interanl diameter (m)*.
opts.L_cc_casing = 609.69e-3;                   % Combustion chamber total casing (pre_cc + cc).
opts.L_pcc = 75e-3;                             % Pre-combustion chamber length.
opts.mass_pcc = 0.5;                            % Pre-combustion chamber mass.
opts.L_cc = 505.8e-3;                           % Combustion chamber total length(m).
opts.T_cc = 3650;                               % Combustion chamber temperature (K).

%% Ox properties.

opts.Molecular_weight_ox = 44.013e-3;           % Molecular weight N2O (kg/mol).
% opts.r_ox = opts.R/opts.Molecular_weight_ox;
opts.gamma_ox = 1.31;                           % Adiabatic index coefficient N2O.
opts.visc_nox = 2.98e-5;                        % Pa.s
opts.calorific_capacity_nox = 2269.5;           % J/kg
opts.thermal_conductivity_nox = 103e-3;         % W/m.K


%% Fuel properties.

opts.L_fuel = 33e-2;            % Fuel length (m).
opts.fuel_mass_init = 3.1;      % Initial fuel mass (kg).
opts.rho_fuel = 900;            % Density of fuel (kg/m^3).
opts.r_fuel_init = 5e-2 / 2;    % Fuel port diameter at ignition.
% opts.r_fuel_init = sqrt(opts.D_cc_int^2/4-opts.fuel_mass_init/(opts.rho_fuel*opts.L_fuel*pi));

opts.reg_a = 15.5e-5;           % Fuel regression parameters of r_dot = a*Gox^n.
opts.reg_n = 0.5;
opts.fuel_margin_mass = 1.2;    % Mass of fuel that is for margin (kg).
opts.fuel_margin_radius = sqrt(opts.D_cc_int^2 / 4 - opts.fuel_margin_mass / (opts.rho_fuel * opts.L_fuel * pi));

opts.CombustionChamberSinusShapeAmplitude = 1/8 ;                                                           % Proportion of initial port radius.
Sin_amp = opts.CombustionChamberSinusShapeAmplitude; 
R = opts.r_fuel_init;
dc = @(theta) sqrt((0.94 * R + R * Sin_amp * sin(8 * theta)).^2 + (R * Sin_amp * 8 * cos(8 * theta)).^2);   % Combustion diameter taking into account sinus shape.
opts.CombustionChamberInitialPerimeter = integral(dc,0,2*pi);                                               % Perimeter taking into account sinus shape.


%% Air properties sea level at 0°.

opts.rho_air_SL = 1.292;                    % Air density (kg/m^3).
opts.visc_dyn_air_SL = 1.729e-5;            % Air dynamic viscosity (kg/m.s).
opts.cp_air_SL = 1006;                      % Specific heat of air (J/kg.K).
opts.air_thermal_conductivity = 0.02364;    % Thermal conductivity air (W/m.K).

%% Combustion properties.

opts.gamma_combustion_products = 1.18;                  % Heat capacity ratio.
opts.molecular_weight_combustion_products = 29e-3;      % Molecular weight of products (kg/mol).
opts.T_cc = 3700;                                       % Combustion temperature (K).


%% Nozzle properties.

opts.D_throat = 38.4e-3;
% opts.A_throat_init = pi*(opts.D_throat)^2/4;  % Nozzle throat area (m^2).
Ae_At = 4.75;
opts.D_exit = sqrt(Ae_At) * opts.D_throat;
opts.A_exit = pi * (opts.D_exit)^2 / 4;         % Nozzle exit area (m^2).

opts.beta_nozzle = 80;                          % Nozzle inlet angle (in °).
opts.alpha_nozzle = 10;                         % Nozzle exit angle (in °).
opts.L_nozzle = 154.55e-3;                      % Nozzle length (m).

%% Tank properties.

opts.aluminium_thermal_conductivity = 236;      % Wm-1K-1 at 0 degree celcius.
opts.rho_alu = 2700;                            % Density aluminium (kg/m^3).
opts.alu_thermal_capacity = 897;                % J/K/kg
opts.aluminium_emissivity_painted = 0.8;        % Emissivity of painted tank.
opts.aluminium_emissivity = 0.3;                % Emissivity of plain aluminium.
opts.aluminium_absorbitivity = 0.4;             % Absorptivity of plain aluminium.

%% Setup the import options.

import_options_N2O = delimitedTextImportOptions("NumVariables", 8);
import_options_c_star = delimitedTextImportOptions("NumVariables", 2);
import_options_CO2 = delimitedTextImportOptions("NumVariables", 8);     % Modified by Benjamin Verbeek 2021-05-11 20:00 CEST "Added CO2 data identically to N2O".

% Specify range and delimiter.
import_options_N2O.DataLines = [8, 602];
import_options_N2O.Delimiter = ";";

import_options_c_star.DataLines = [2,23];
import_options_c_star.Delimiter = ";";

import_options_CO2.DataLines = [8, 602];
import_options_CO2.Delimiter = ";";

% Specify column names and types.
import_options_N2O.VariableNames = ["TemperatureK", "Pressurebar", "Liquiddensitykgm", "Gasdensitykgm", "LiquidIntEnergy", "VaporIntEnergy", "LiquidEnthalpy", "VaporEnthalpy"];
import_options_N2O.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];
import_options_N2O.ExtraColumnsRule = "ignore";
import_options_N2O.EmptyLineRule = "read";

import_options_c_star.VariableNames = ["OF", "CStarms"];
import_options_c_star.VariableTypes = ["double", "double"];
import_options_c_star.ExtraColumnsRule = "ignore";
import_options_c_star.EmptyLineRule = "read";

import_options_CO2.VariableNames = ["TemperatureK", "Pressurebar", "Liquiddensitykgm", "Gasdensitykgm", "LiquidIntEnergy", "VaporIntEnergy", "LiquidEnthalpy", "VaporEnthalpy"];
import_options_CO2.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];
import_options_CO2.ExtraColumnsRule = "ignore";
import_options_CO2.EmptyLineRule = "read";

% Import the data.
NO2 = readtable("./datasets/nitrous-oxide_LVsaturation.csv", import_options_N2O);
C_star = readtable("./datasets/characteristic_velocity.csv", import_options_c_star);
CO2 = readtable("./datasets/carbon-dioxide_LVsaturation.csv", import_options_CO2);

%% Clear temporary variables.

clear import_options_N2O
clear import_options_c_star
clear import_options_CO2

%% Spline fitting.

Temperature_set = NO2.TemperatureK;         % Get temperature range.
NO2_Psat_set = NO2.Pressurebar;             % Get saturation pressure for the temperatures above.
NO2_Rhol_set = NO2.Liquiddensitykgm;        % Get liquid density for the temperatures above.
NO2_Rhog_set = NO2.Gasdensitykgm;           % Get gas density for the temperatures above.
NO2_Ul_set = NO2.LiquidIntEnergy;           % Get liquid internal energy for the temperatures above.
NO2_Ug_set = NO2.VaporIntEnergy;            % Get gas internal energy for the temperatures above.

% Put all variables in a map under the shorthand names.
names = {'T', 'Psat', 'Rhol', 'Rhog', 'Ul', 'Ug'};
vars = [Temperature_set, NO2_Psat_set, NO2_Rhol_set, NO2_Rhog_set, NO2_Ul_set, NO2_Ug_set];
NO2_vars = containers.Map;
for i = 1:length(names)
    NO2_vars(string(names(i))) = vars(:, i);
end

% Specify the x and y variables that are required, and a name for each pair.
xs = {'T', 'T', 'T', 'T', 'T', 'Psat', 'Psat', 'Psat'};
ys = {'Psat', 'Rhol', 'Rhog', 'Ul', 'Ug', 'Rhog', 'Ul', 'Ug'};
names = {'Psat', 'RhoL_T', 'RhoG_T', 'UL_T', 'UG_T', 'RhoG_P', 'UL_P', 'UG_P'};

% Fit a spline to each (x, y) pair.
for i = 1:length(names)
    x = NO2_vars(string(xs(i)));
    y = NO2_vars(string(ys(i)));
    name = string(names(i)) + '_NO2_spline';
    spln = csaps(x, y);                         % Fit a cubic smoothing spline to the data.
    opts.(name) = fnxtr(spln, 2);               % Extrapolate with a quadratic polynomial to avoid wonkiness at the boundaries.
end

% Plots for debugging.
debug_plot = false;
if debug_plot
    rows = ceil(sqrt(length(names)));
    tiledlayout(rows, rows)
    for i = 1:length(names)
        x = NO2_vars(string(xs(i)));
        y = NO2_vars(string(ys(i)));
        name = string(names(i)) + '_NO2_spline';
        nexttile
        hold on
        xlim([0.99 * min(x) 1.01 * max(x)])
        scatter(x, y, '.')
        fnplt(opts.(name))
        hold off
    end
end

% Delete temporary variables.
clear debug_plot i name names NO2_vars spln vars x xs y ys;

opts.OF_set = C_star.OF;                                % OF ratio range.
opts.C_star_set = C_star.CStarms;                       % Characteristic velocity C_Star.
% opts.C_Star_polynom=polyfit(OF_set,C_star_set,5);     % Interpolation degree 3.

Temperature_set = CO2.TemperatureK;         % Get temperature range.
CO2_Psat_set = CO2.Pressurebar;             % Get saturation pressure for the temperatures above.
CO2_Rhol_set = CO2.Liquiddensitykgm;        % Get liquid density for the temperatures above.
CO2_Rhog_set = CO2.Gasdensitykgm;           % Get gas density for the temperatures above.
CO2_Ul_set = CO2.LiquidIntEnergy;           % Get liquid internal energy for the temperatures above.
CO2_Ug_set = CO2.VaporIntEnergy;            % Get gas internal energy for the temperatures above.


%% Storage tank geometry.

% TODO: Make sure that this data is only used for the tank filling
%       simulation and remove/move it.

opts.D_ext_storage = 230e-3;                                        % Storage tank external diameter (m).
opts.V_storage = 50e-3;                                             % Storage tank volume (m^3).
opts.D_int_storage = opts.D_ext_storage - 2 * opts.e_tank;
opts.L_storage = opts.V_storage / (pi * (opts.D_int_storage)^2 / 4);

%% Filling properties.

% TODO: Make sure that this data is only used for the tank filling
%       simulation and remove/move it.

opts.d_filling_inlet = 4.7e-3;      %2.5e-3;    % m 
opts.d_filling_outlet = 0.9e-3;                 % m

opts.S_inlet = pi * (opts.d_filling_inlet)^2 / 4;
opts.S_outlet = pi * (opts.d_filling_outlet)^2 / 4;

opts.P_storage_tank_init = fnval(opts.Psat_NO2_spline, opts.T_ext) * 10^5;     % NOTE: not changed /Benjamin.
opts.cd_inlet = 0.85;
opts.cd_outlet = 0.95;
% opts.r_ox = py.CoolProp.CoolProp.PropsSI('P','T',opts.T_ext,'Q', 1,'NitrousOxide') / py.CoolProp.CoolProp.PropsSI('D','T',opts.T_ext,'Q', 1,'NitrousOxide') / opts.T_ext;
% opts.r_ox = 180.7175;
