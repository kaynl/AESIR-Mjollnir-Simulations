function m_dot_th = Mass_flow_throat(P_cc,OF,At)
%MASS_FLOW_TH calculates the mass flow of that goes through the throat

global opts
c_star = interp1q(opts.OF_set,opts.C_star_set,OF);

m_dot_th = P_cc*At/c_star;

end

