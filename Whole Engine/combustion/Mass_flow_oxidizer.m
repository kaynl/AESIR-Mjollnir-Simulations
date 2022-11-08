function mf_ox = Mass_flow_oxidizer(T_tank, P_tank, P_cc)
%MASS_FLOW_FUEL calculates the mass flow of the fuel
global opts

[P_cc_range, mf_crit] = critical_mf_MOODY_NO2(P_tank, T_tank, opts.Cd);

mf_ox = interp1(P_cc_range, mf_crit, P_cc);

end

