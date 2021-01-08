function [critical_mf,critical_P_cc] = critical_mf_DYER(P_tank,T_tank)
%CRITICAL_MF_DYER Summary of this function goes here
%   Detailed explanation goes here

    global opts
    Cd = opts.Cd;                     %Discharge coefficient
    D_inj = 2*opts.r_inj;             %Injector diameter
    n_inj = opts.n_inj;               %Number of injector holes
    Ai=n_inj*pi*D_inj^2/4;
    
    P_cc = 1e5:1e5:P_tank;
    mf_SPI = zeros(1,length(P_cc));
    mf_HEM = zeros(1,length(P_cc));
    kapa = ones(1,length(P_cc));
    
    h1 = py.CoolProp.CoolProp.PropsSI('H','P',P_tank,'T', round(T_tank,2),'NitrousOxide');   %Enthalpie massic (J/kg)
    s1 = py.CoolProp.CoolProp.PropsSI('S','P',P_tank,'T', round(T_tank,2),'NitrousOxide');   %Enthalpie massic (J/K.kg)
    
    P_sat = polyval(opts.Psat_NO2_polynom,T_tank)*10^6;
    rho_Ox_1 = polyval(opts.RhoL_T_NO2_polynom,T_tank);           %Density of Oxidizer (kg/m^3)
    rho_Ox_2 = polyval(opts.RhoL_Psat_NO2_polynom,P_cc/10^6);     %Density of Oxidizer (kg/m^3)
    
    
    for i=1:length(P_cc)
        if P_sat>P_cc(i)
            kapa(i) = sqrt((P_tank-P_cc)/(P_sat-P_cc));
        else
            kapa(i) = 1;
        end
        h2 = py.CoolProp.CoolProp.PropsSI('H','P',P_cc(i),'S', s1,'NitrousOxide');   %Enthalpie massic (J/kg)
        
        mf_SPI(i) = Cd*Ai*sqrt(2*rho_Ox_1*(P_tank-P_cc(i)));
        mf_HEM(i) = Cd*Ai*rho_Ox_2(i)*sqrt(2*(h1-h2));
        
    end
    
    
    mf_ox = (kapa.*mf_SPI+mf_HEM./(1+kapa));
    
    figure(1)
    plot(P_tank-P_cc,mf_ox,P_tank-P_cc,mf_SPI,P_tank-P_cc,mf_HEM)
    legend("Total","SPI","HEM")
    xlabel("Pressure drop (Pa)")
    ylabel("Mass flow (kg/s)")
    title("DYER Mass flow (P1=59 bars, T1 = 287K)")
    
    [critical_mf,index_critical_P_cc] = max(mf_ox);
    
    critical_P_cc = P_cc(index_critical_P_cc);
    
end

