function x_vapor = x_vapor(U_tot,m_tot,T)
%This function computes the value of the vapor factor tanks to the
%temperatre inside the tank

    %Coolprop u(T) liquid at saturation (J/kg)
    u_liq = py.CoolProp.CoolProp.PropsSI('U','T',T,'Q', 0,'NitrousOxide');
    
    %Coolprop u(T) vapor at saturation (J/kg)
    u_vap = py.CoolProp.CoolProp.PropsSI('U','T',T,'Q', 1,'NitrousOxide');
    
    
    x_vapor=((U_tot./m_tot)-u_liq)./(u_vap-u_liq);

end

