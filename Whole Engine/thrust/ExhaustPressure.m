function Pe = ExhaustPressure(Pcc,P_ext,Me,opts)
%EXHAUSTPRESSURE Calculates the exhaust pressure

gam = opts.gamma_combustion_products;

if Pcc == P_ext
    Pe=P_ext;
else
    Pe = Pcc.*(1+(gam-1)*Me^2/2)^(gam/(1-gam));
end

end

