function plot_metric_distribution(metric_values, dynamic_threshold, metric_name, clusterName, plotDir, suffix)
    valid_mask = metric_values > 0;
    surviving_values = metric_values(valid_mask);
    
    f = figure('Visible', 'off');
    if ~isempty(surviving_values)
        log_vals = log10(surviving_values);
        histogram(log_vals, 50, 'Normalization', 'pdf');
        hold on; 
        
        % FIX: Ensure threshold uses log10 to match the data
        log_thresh = log10(dynamic_threshold);
        xl = xline(log_thresh, '--r', 'Threshold', 'LineWidth', 2, 'LabelOrientation', 'aligned');
        xl.FontSize = 10;
        
        title(sprintf('Log10 Density (%s): %s', clusterName, suffix), 'Interpreter', 'none');
        xlabel(['log10(', metric_name, ')']);
        ylabel('Probability Density');
        hold off;
    else
        title(sprintf('No %s > 0 found (%s)', metric_name, clusterName), 'Interpreter', 'none');
        text(0.5, 0.5, 'Matrix is empty', 'HorizontalAlignment', 'center');
    end
    
    saveas(f, fullfile(plotDir, sprintf('%s_%s_%s.png', clusterName, metric_name, suffix)));
    close(f);
end