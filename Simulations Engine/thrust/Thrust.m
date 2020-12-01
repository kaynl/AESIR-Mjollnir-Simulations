function [Tr] = Thrust(mf_throat,Ve,P_ext,Pe,opts)
%THRUST of the rocket

Ae = opts.A_exit;

Tr = mf_throat*Ve+(Pe-P_ext)*Ae;

end

