function oim_phase_plotter(confs, ts)
    sz = size(confs); sz = sz(1:end-1);
    phases_by_oscillator = reshape(confs,[prod(sz) numel(confs)/prod(sz)]);
    
    figure;
    hold on;
    for i = 1:size(phases_by_oscillator, 1)
        plot(log(ts), phases_by_oscillator(i, :));
    end
    
    % Add a horizontal black thick dashed line at pi/2
    yline(pi/2, '--k', 'LineWidth', 2); % '--k' specifies a black dashed line, 'LineWidth', 2 specifies the thickness
    
    hold off;
    xlabel('log(time)');
    ylabel('Phases');
    title('Oscillator Phases Over Time');

    % Save the figure
    saveas(gcf, 'oim_phase_plot.png');

    % Close the figure
    close(gcf);
end
