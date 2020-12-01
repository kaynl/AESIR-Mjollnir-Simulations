function f_Mach = ExhaustMach_fct(Mach,At,opts)
%EXHAUSTMACH_FCT fct to find the zero in order to find the exhaust mach
%number

Ae = opts.A_exit;
gam = opts.gamma_combustion_products;

f_Mach = 1./Mach.*(2/(gam+1).*(1+(gam-1).*Mach.^2./2)).^((gam+1)/(2*(gam-1)))-Ae/At;


end

