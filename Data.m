%%% This file gathers all engine data needed for the simulation

% TODO: Split data into parameters and "constants" based on how we need
%       to use them.
% TODO: Create a good way to run a specific simulation with a certain set
%       of parameter choices.

global opts

opts.flight_state = 1;          %0 if test fire on the grouond and 1 if launch
opts.filling_ratio = 0.95;         %Tank filling ratio
opts.launch_angle = 87;          %Self explanatory (in °)

opts.drag_coefficient = 0.5;     
opts.combustion_efficiency = 0.9;
opts.T_ext = 293;               %(K) Exterior temperature (20°C)

%% Physical Constants

opts.g=9.81;                                    %Gravitational Constant (m.s-2)
opts.R = 8.314;                                 %Universal Gas Constant (J⋅K−1⋅mol−1) 

opts.P_atm_sl = 101325;                         %Atmospheric Pressure (N/m2)
opts.stephan_cst = 5.67e-8;                     %Stephan-Boltzman Constant (W/m2/K4)
opts.eber_parameter = 0.89;                     %Eber Parameter for vertex angle between 20-50 degrees
opts.Molecular_weight_air = 28.9647e-3;         %Molecular Weight of Air (kg/mol)
opts.r_air = opts.R/opts.Molecular_weight_air;


%% Requirements

opts.design_altitude = 14000;          %Designed Altitude to reach (m)
opts.required_altitude = 12000;        %Mission requirements (m)


%% Mass
opts.parachute_mass = 10;
opts.electronics_mass = 2.3;
opts.bodyTube_mass = 7;
opts.payload_mass = 2;

opts.propulsionSystem = 24.504;

opts.dry_mass = opts.parachute_mass + opts.electronics_mass + opts.bodyTube_mass + opts.payload_mass + opts.propulsionSystem;

opts.m_ox_init=24.5;              %Oxidizer mass (kg)
opts.m_fuel_init=3.1;             %Fuel Mass (kg)
opts.rho_ox = 785;                %Oxidizer density (kg/m^3)

%% Tank Geometry

opts.D_ext_tank = 16e-2;%10e-2;%    %Tank external diameter (m)
opts.e_tank = 3.5e-3;          %Tank thickness
opts.D_int_tank = opts.D_ext_tank-2*opts.e_tank;%9.42e-2;%    %Tank internal diameter (m)
opts.L_tank = 1.83;%0.73;%      %Tank Length (m)
opts.V_tank = pi*(opts.D_int_tank)^2/4*opts.L_tank;%33.1e-3;  %      %Tank Volume (m^3) (present in Tank_Temperature_finder_fct)
opts.surface = pi*(opts.D_ext_tank)^2/4; %Rocket Surface

%% Kastrullen
opts.L_kastrullen = 35e-2;  %length of Kastrullen

%% Injector Geometry
              
opts.r_inj=1.2e-3/2;           %injector radius (m)
opts.Cd = 0.89;                %Discharge coefficient
opts.L_inj = 15e-3;            %Injector Plate thickness (m)
opts.n_inj = 34;               %Number of injectors

opts.r_inj_plate = 30e-3;      %m
opts.mass_inj = 0.271;         %kg
opts.e_inj = 0.013;             %m

%% Combustion Chamber Geometry

opts.D_cc_ext = 15.2e-2;                      %Combustion Chamber external diameter (m)
opts.e_cc = 4e-3;
opts.D_cc_int = opts.D_cc_ext-2*opts.e_cc;     %Combustion Chamber interanl diameter (m)*
opts.L_cc_casing = 609.69e-3;               %Combustion Chamber Total Casing (pre_cc + cc)
opts.L_pcc = 75e-3;                      %Pre-combustion chamber length
opts.mass_pcc = 0.5;                     %Pre-combustion chamber mass
opts.L_cc = 505.8e-3;                          %Combustion Chamber Total length(m)
opts.T_cc = 3650;                           %Combustion Chamber temperature (K)

%% Ox Properties

opts.Molecular_weight_ox = 44.013e-3;  %molecular weight N2O (kg/mol)
% opts.r_ox = opts.R/opts.Molecular_weight_ox;
opts.gamma_ox = 1.31;       %Adiabatic Index Coefficient N2O
opts.visc_nox = 2.98e-5;   %Pa.s
opts.calorific_capacity_nox = 2269.5;       %J/kg
opts.thermal_conductivity_nox = 103e-3;             %W/m.K


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
opts.T_cc = 3700;                                       %Combustion temperature (K)


%% Nozzle Properties

opts.D_throat =38.4e-3;
%opts.A_throat_init = pi*(opts.D_throat)^2/4;      %Nozzle Throat Area (m^2)
Ae_At = 4.75;
opts.D_exit = sqrt(Ae_At)*opts.D_throat;
opts.A_exit = pi*(opts.D_exit)^2/4;          %Nozzle Exit Area (m^2)

opts.beta_nozzle = 80;                       %Nozzle Inlet Angle (in °)
opts.alpha_nozzle = 10;                      %Nozzle Exit Angle (in °)
opts.L_nozzle = 154.55e-3;                   %Nozzle Length (m)

%% Tank Properties

opts.aluminium_thermal_conductivity=236;    %Wm-1K-1 at 0 degree celcius
opts.rho_alu = 2700;                        %Density Aluminium (kg/m^3)
opts.alu_thermal_capacity = 897;            %J/K/kg
opts.aluminium_emissivity_painted = 0.8;    %Emissivity of painted tank
opts.aluminium_emissivity = 0.3;           %Emissivity of plain aluminium
opts.aluminium_absorbitivity = 0.4;         %Absorptivity of plain aluminium

%% Setup the Import Options
import_options_N2O = delimitedTextImportOptions("NumVariables", 8);
import_options_c_star = delimitedTextImportOptions("NumVariables", 2);
import_options_CO2 = delimitedTextImportOptions("NumVariables", 8); % Modified by Benjamin Verbeek 2021-05-11 20:00 CEST "Added CO2 data identically to N2O"

% Specify range and delimiter
import_options_N2O.DataLines = [8, 602];
import_options_N2O.Delimiter = ";";

import_options_c_star.DataLines = [2,23];
import_options_c_star.Delimiter = ";";

import_options_CO2.DataLines = [8, 602];
import_options_CO2.Delimiter = ";";

% Specify column names and types
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

% Import the data
NO2 = readtable("./datasets/nitrous-oxide_LVsaturation.csv", import_options_N2O);
C_star = readtable("./datasets/characteristic_velocity.csv", import_options_c_star);
CO2 = readtable("./datasets/carbon-dioxide_LVsaturation.csv", import_options_CO2);

%% Clear temporary variables
clear import_options_N2O
clear import_options_c_star
clear import_options_CO2

%% INTERPOLATION

Temperature_set=NO2.TemperatureK;                 %Getting temperature range
NO2_Psat_set=NO2.Pressurebar;                     %Getting saturation pressure for the temperatures above
NO2_Rhol_set=NO2.Liquiddensitykgm;                %Getting liquid density for the temperatures above
NO2_Rhog_set=NO2.Gasdensitykgm;                   %Getting gaz density for the temperatures above
NO2_Ul_set=NO2.LiquidIntEnergy;                   %Getting liquid internal energy for the temperatures above
NO2_Ug_set=NO2.VaporIntEnergy;                    %Getting gaz internal Energy for the temperatures above

% TODO: Make sure that degree 3 gives a good interpolation
%       A: doesn't really seem so.
% TODO: Investigate how and where the interpolation results are used
% TODO: Update to a better interpolation method where necessary (eg.
%       linear)

opts.Psat_NO2_polynom=polyfit(Temperature_set,NO2_Psat_set,3);  %interpolation polynomial of degree 3
opts.RhoL_Psat_NO2_polynom=polyfit(NO2_Psat_set,NO2_Rhol_set,3);    %interpolation polynomial of degree 3

opts.RhoL_T_NO2_polynom=polyfit(Temperature_set,NO2_Rhol_set,3);    %interpolation polynomial of degree 3
opts.RhoG_T_NO2_polynom=polyfit(Temperature_set,NO2_Rhog_set,3);    %interpolation polynomial of degree 3
opts.RhoL_P_NO2_polynom=polyfit(NO2_Psat_set,NO2_Rhol_set,3);    %interpolation polynomial of degree 3
opts.RhoG_P_NO2_polynom=polyfit(NO2_Psat_set,NO2_Rhog_set,3);    %interpolation polynomial of degree 3

opts.UL_T_NO2_polynom=polyfit(Temperature_set,NO2_Ul_set,3);    %interpolation polynomial of degree 3
opts.UG_T_NO2_polynom=polyfit(Temperature_set,NO2_Ug_set,3);    %interpolation polynomial of degree 3
opts.UL_P_NO2_polynom=polyfit(NO2_Psat_set,NO2_Ul_set,3);    %interpolation polynomial of degree 3
opts.UG_P_NO2_polynom=polyfit(NO2_Psat_set,NO2_Ug_set,3);    %interpolation polynomial of degree 3


opts.OF_set = C_star.OF;                                             %OF ratio range
opts.C_star_set = C_star.CStarms;                                    %characteristic velocity C_Star
% opts.C_Star_polynom=polyfit(OF_set,C_star_set,5);                  %interpolation degree 3



Temperature_set=CO2.TemperatureK;                 %Getting temperature range
CO2_Psat_set=CO2.Pressurebar;                     %Getting saturation pressure for the temperatures above
CO2_Rhol_set=CO2.Liquiddensitykgm;                %Getting liquid density for the temperatures above
CO2_Rhog_set=CO2.Gasdensitykgm;                   %Getting gaz density for the temperatures above
CO2_Ul_set=CO2.LiquidIntEnergy;                   %Getting liquid internal energy for the temperatures above
CO2_Ug_set=CO2.VaporIntEnergy;                    %Getting gaz internal Energy for the temperatures above


opts.Psat_CO2_polynom=polyfit(Temperature_set,CO2_Psat_set,3);  %interpolation polynomial of degree 3
opts.RhoL_Psat_CO2_polynom=polyfit(CO2_Psat_set,CO2_Rhol_set,3);    %interpolation polynomial of degree 3

opts.RhoL_T_CO2_polynom=polyfit(Temperature_set,CO2_Rhol_set,3);    %interpolation polynomial of degree 3
opts.RhoG_T_CO2_polynom=polyfit(Temperature_set,CO2_Rhog_set,3);    %interpolation polynomial of degree 3
opts.RhoL_P_CO2_polynom=polyfit(CO2_Psat_set,CO2_Rhol_set,3);    %interpolation polynomial of degree 3
opts.RhoG_P_CO2_polynom=polyfit(CO2_Psat_set,CO2_Rhog_set,3);    %interpolation polynomial of degree 3

opts.UL_T_CO2_polynom=polyfit(Temperature_set,CO2_Ul_set,3);    %interpolation polynomial of degree 3
opts.UG_T_CO2_polynom=polyfit(Temperature_set,CO2_Ug_set,3);    %interpolation polynomial of degree 3
opts.UL_P_CO2_polynom=polyfit(CO2_Psat_set,CO2_Ul_set,3);    %interpolation polynomial of degree 3
opts.UG_P_CO2_polynom=polyfit(CO2_Psat_set,CO2_Ug_set,3);    %interpolation polynomial of degree 3

%% Storage Tank Geometry

% TODO: Make sure that this data is only used for the tank filling
%       simulation and remove/move it.

opts.D_ext_storage = 230e-3;        %Storage Tank external diameter (m)
opts.V_storage = 50e-3;             %Storage Tank Volume (m^3)
opts.D_int_storage = opts.D_ext_storage - 2*opts.e_tank;
opts.L_storage = opts.V_storage/(pi*(opts.D_int_storage)^2/4);


%% Filling Properties

% TODO: Make sure that this data is only used for the tank filling
%       simulation and remove/move it.

opts.d_filling_inlet = 4.7e-3;%2.5e-3;%m 
opts.d_filling_outlet = 0.9e-3;%m

opts.S_inlet = pi*(opts.d_filling_inlet)^2/4;
opts.S_outlet = pi*(opts.d_filling_outlet)^2/4;


opts.P_storage_tank_init = polyval(opts.Psat_NO2_polynom,opts.T_ext)*10^5;  % NOTE: not changed /Benjamin
opts.cd_inlet = 0.85;
opts.cd_outlet = 0.95;
%opts.r_ox = py.CoolProp.CoolProp.PropsSI('P','T',opts.T_ext,'Q', 1,'NitrousOxide') / py.CoolProp.CoolProp.PropsSI('D','T',opts.T_ext,'Q', 1,'NitrousOxide') / opts.T_ext;
% opts.r_ox = 180.7175;