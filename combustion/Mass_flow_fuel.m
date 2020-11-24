function m_dot_fuel = Mass_flow_fuel(G_ox,r_fuel)
%MASS_FLOW_FUEL calculates the mass flow of the fuel

global opts

rho_fuel = opts.rho_fuel;
L = opts.L_fuel;
a = opts.reg_a;
n = opts.reg_n;
r_dot = a*(G_ox.^n);

Sin_amp = opts.CombustionChamberSinusShapeAmplitude;
C8 = opts.CombustionChamberInitialPerimeter;
Rinit=opts.r_fuel_init;
coeff = 1+(C8/(2*pi*Rinit) - 1)*exp(sqrt(6/Sin_amp)*(Rinit-r_fuel)*2/Rinit);

m_dot_fuel = rho_fuel*2*pi*r_fuel*r_dot*L*coeff;

end

