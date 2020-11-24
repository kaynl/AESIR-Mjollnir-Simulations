function [V_e] = ExhaustSpeed(T_cc,P_e,P_cc,opts)
%EJECTION_SPEED Gaz ejection speed calculator based on Laval Nozzle formula

R=opts.R;
M=opts.Molecular_weigth_combustion_products;
gamma = opts.gamma_combustion_products;

V_e = sqrt((T_cc.*R/M)*2*gamma./(gamma-1).*(1-(P_e./P_cc).^((gamma-1)/(gamma))));

end

