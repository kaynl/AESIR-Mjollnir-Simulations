function ht_sensors_compilation(ht_data, test_name, t0, t1)
    sample_points = linspace(t0, t1, 1000);
    
    rows = 3;
    columns = 1;
    pressure_range = [0 50];
    figure();
    sgtitle(test_name)
    
    % feed pressures
    subplot(rows, columns, 1);
    plot(sample_points, ht_data.PTRAN_1_I(sample_points));
    hold on
    plot(sample_points, ht_data.PTRAN_2_I(sample_points));
    plot(sample_points, ht_data.PTRAN_3_I(sample_points));
    lgd = legend("1: Top Pressure", "2: Bottom Pressure", "3: Injector Pressure");
    ylim(pressure_range);
    grid on
    title("Feed Pressure")
    xlabel("Time [s]")
    ylabel("Pressure [bar]")
    hold off
    
    % combustion pressure
    subplot(rows, columns, 2);
    plot(sample_points, ht_data.PTRAN_4_I(sample_points));
    hold on
    plot(sample_points, ht_data.PTRAN_1_I(sample_points), '--');
    lgd = legend("4: Combustion Pressure", "1: Top Pressure (REF)");
    ylim(pressure_range);
    grid on
    title("Combustion Pressure")
    xlabel("Time [s]")
    ylabel("Pressure [bar]")
    hold off
    
    % thrust
    subplot(rows, columns, 3);
    plot(sample_points, ht_data.THRUST_I(sample_points));
    grid on
    title("Thrust")
    xlabel("Time [s]")
    ylabel("Force [kN]")
    
    set(gcf, 'Position', [0 0 1000 1000]);
end