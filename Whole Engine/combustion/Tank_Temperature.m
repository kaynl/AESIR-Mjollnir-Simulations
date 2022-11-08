function T_tank = Tank_Temperature(U_total,m_total)
%TANK_TEMPERATURE : finding the zero of the thermodynamic equilibrium function to get internal temp of the tank
%disp(U_total)
%disp(m_total)
T_tank = fzero(@(T_unknown) Tank_Temperature_finder_fct(U_total,m_total,T_unknown), [183, 309.51]);

end

