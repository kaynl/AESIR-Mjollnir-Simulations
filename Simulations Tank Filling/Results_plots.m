T_ext=opts.T_ext;
T_tank=T_ext;
rho_liq=py.CoolProp.CoolProp.PropsSI('D','T',T_tank,'Q', 0,'NitrousOxide');

plot(t/60,100*m_liq/opts.V_tank/rho_liq,t/60,100-100*m_liq/opts.V_tank/rho_liq)
ylabel('Tank filling in %')
xlabel('Time in min')
title('Tank filling versus time')
legend('V_l_i_q','V_v_a_p (*)')%(*) atm simple version of V_vap