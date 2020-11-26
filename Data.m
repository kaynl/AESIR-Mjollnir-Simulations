%%% This file gathers all engine data needed for the simulation

global opts

opts.flight_state = 1;          %0 if test fire on the grouond and 1 if launch
opts.filling_ratio = 1;         %Tank filling ratio
opts.launch_angle = 85;          %Self explanatory (in °)

opts.drag_coefficient = 0.5;     
opts.combustion_efficiency = 0.9;

%% Requirements

opts.design_altitude = 14000;          %Designed Altitude to reach (m)
opts.required_altitude = 12000;        %Mission requirements (m)


%% Mass
opts.parachute_mass = 10;
opts.electronics_mass = 2.3;
opts.bodyTube_mass = 7;
opts.payload_mass = 2;
opts.dry_mass = opts.parachute_mass + opts.electronics_mass + opts.bodyTube_mass + opts.payload_mass;

opts.m_ox_init=24.5;              %Oxidizer mass (kg)
opts.m_fuel_init=3.1;             %Fuel Mass (kg)
opts.rho_ox = 785;                %Oxidizer density (kg/m^3)

%% Tank Geometry

opts.D_tank_ext = 16e-2;       %Tank external diameter (m)
opts.D_tank_int = 15.42e-2;    %Tank internal diameter (m)
opts.e_tank = 2.9e-3;          %Tank thickness
opts.L_tank = 1.73;            %Tank Length (m)
opts.V_tank = 31.2e-3;         %Tank Volume (m^3) (present in Tank_Temperature_finder_fct)
opts.surface = pi*(opts.D_tank_ext)^2/4; %Rocket Surface

%% Kastrullen
opts.L_kastrullen = 35e-2;  %length of Kastrullen

%% Injector Geometry

opts.n_inj=38;                 %Number of injectors
opts.r_inj=1.5e-3/2;             %injector radius (m)
opts.Cd = 0.83;                %Discharge coefficient ()

%% Combustion Chamber Geometry

opts.D_cc_ext = 15e-2;                      %Combustion Chamber external diameter (m)
opts.D_cc_int = opts.D_cc_ext-2*2.1e-3;     %Combustion Chamber interanl diameter (m)
opts.L_pcc = 103.9e-3;                      %Pre-combustion chamber length
opts.L_cc = 50e-2;                          %Combustion Chamber Total length(m)
opts.T_cc = 3500;                           %Combustion Chamber temperature (K)


%% Fuel Properties

opts.L_fuel = 33e-2;            %Fuel length (m)
opts.fuel_mass_init = 3.1;      %Initial fuel mass (kg)
opts.rho_fuel=900;              %Density of fuel (kg/m^3)
opts.r_fuel_init = 5e-2/2;      %fuel port diameter at ignition
% opts.r_fuel_init = sqrt(opts.D_cc_int^2/4-opts.fuel_mass_init/(opts.rho_fuel*opts.L_fuel*pi));

opts.reg_a = 15.5e-5;          %Fuel Regression parameters of r_dot = a*Gox^n
opts.reg_n = 0.5;
opts.fuel_margin_mass = 1.2;    %mass of fuel that is for margin (kg)
opts.fuel_margin_radius = sqrt(opts.D_cc_int^2/4-opts.fuel_margin_mass/(opts.rho_fuel*opts.L_fuel*pi));

opts.CombustionChamberSinusShapeAmplitude = 1/8 ; % Proportion of initial port radius
Sin_amp = opts.CombustionChamberSinusShapeAmplitude; 
R=opts.r_fuel_init;
dc = @(theta) sqrt((0.94*R+R*Sin_amp*sin(8*theta)).^2+(R*Sin_amp*8*cos(8*theta)).^2); %combustion diameter taking into account sinus shape
opts.CombustionChamberInitialPerimeter = integral(dc,0,2*pi); %perimeter taking into account sinus shape


%% Air Properties sea level at 0°

opts.rho_air_SL = 1.292;                    %Air density (kg/m^3)
opts.visc_dyn_air_SL = 1.729e-5;            %Air dynamic viscosity (kg/m.s)
opts.cp_air_SL = 1006;                      %Specific Heat of Air (J/kg.K)
opts.air_thermal_conductivity = 0.02364;    %Thermal Conductivity Air (W/m.K)

%% Combustion Properties

opts.gamma_combustion_products = 1.18;                   %Adiabatic Index Coefficient
opts.Molecular_weigth_combustion_products = 29e-3;      %Molecular Weigth of products (kg/mol)
opts.T_cc = 3500;                   %Combustion temperature (K)


%% Nozzle Properties

opts.D_throat = 39.37e-3;
% opts.A_throat_init = pi*(opts.D_throat)^2/4;      %Nozzle Throat Area (m^2)
opts.D_exit = 89.55e-3;
opts.A_exit = pi*(opts.D_exit)^2/4;          %Nozzle Exit Area (m^2)

opts.beta_nozzle = 80;                       %Nozzle Inlet Angle (in °)
opts.alpha_nozzle = 10;                      %Nozzle Exit Angle (in °)
opts.L_nozzle = 154.55e-3;                   %Nozzle Length (m)

%% Tank Properties

opts.aluminium_thermal_conductivity=236;    %Wm-1K-1 at 0 degree celcius
opts.rho_alu = 2700;                        %Density Aluminium (kg/m^3)
opts.alu_thermal_capacity = 897;            %J/K/kg
opts.aluminium_emissivity = 0.8;            %Emissivity of painted tank

%% Physical Constants

opts.g=9.81;                   %Gravitational Constant (m.s-2)
opts.R = 8.314;                %Universal Gas Constant (J⋅K−1⋅mol−1)    
opts.P_atm_sl = 101325;        %Atmospheric Pressure (N/m2)
opts.stephan_cst = 5.67e-8;    %Stephan-Boltzman Constant (W/m2/K4)
opts.eber_parameter = 0.89;    %Eber Parameter for vertex angle between 20-50 degrees

%% Setup the Import Options
import_options_N2O = delimitedTextImportOptions("NumVariables", 4);
import_options_c_star = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
import_options_N2O.DataLines = [2, 257];
import_options_N2O.Delimiter = ";";

import_options_c_star.DataLines = [2,23];
import_options_c_star.Delimiter = ";";

% Specify column names and types
import_options_N2O.VariableNames = ["TemperatureK", "Pressurebar", "Liquiddensitykgm", "Gasdensitykgm"];
import_options_N2O.VariableTypes = ["double", "double", "double", "double"];
import_options_N2O.ExtraColumnsRule = "ignore";
import_options_N2O.EmptyLineRule = "read";

import_options_c_star.VariableNames = ["OF", "CStarms"];
import_options_c_star.VariableTypes = ["double", "double"];
import_options_c_star.ExtraColumnsRule = "ignore";
import_options_c_star.EmptyLineRule = "read";

% Import the data
NO2 = readtable("./datasets/nitrous-oxide_LVsaturation.csv", import_options_N2O);
C_star = readtable("./datasets/characteristic_velocity.csv", import_options_c_star);

%% Clear temporary variables
clear import_options_N2O
clear import_options_c_star

%% INTERPOLATION

Temperature_set=NO2.TemperatureK;                 %Getting temperature range
NO2_Psat_set=NO2.Pressurebar;                     %Getting saturation pressure for the temperatures above
NO2_Rho_set=NO2.Liquiddensitykgm;                 %Getting density for the temperatures above

opts.Psat_NO2_polynom=polyfit(Temperature_set,NO2_Psat_set,3);  %interpolation polynomial of degree 3
opts.Rho_T_NO2_polynom=polyfit(Temperature_set,NO2_Rho_set,3);    %interpolation polynomial of degree 3
opts.Rho_Psat_NO2_polynom=polyfit(NO2_Psat_set,NO2_Rho_set,3);    %interpolation polynomial of degree 3



opts.OF_set = C_star.OF;                                             %OF ratio range
opts.C_star_set = C_star.CStarms;                                    %characteristic velocity C_Star
% opts.C_Star_polynom=polyfit(OF_set,C_star_set,5);               %interpolation degree 3
