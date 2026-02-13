function plot_efm_length_distribution(E, saveDir, clusterName)

    
    % Calculate lengths
    efm_lengths = sum(abs(E) > 1e-9, 1);
    
    % Create Plot
    fig = figure('Name', 'EFM Length Analysis', 'Color', 'w', 'Visible', 'off');
    
    % Histogram
    h = histogram(efm_lengths);
    h.FaceColor = [0.2 0.6 0.5];
    h.EdgeColor = 'w';
    h.FaceAlpha = 0.8;
    
    % Formatting
    grid on;
    xlabel('EFM Length (Number of Active Reactions)', 'FontSize', 11);
    ylabel('Count (Number of EFMs)', 'FontSize', 11);
    
    % Statistics
    avg_len = mean(efm_lengths);
    med_len = median(efm_lengths);
    min_len = min(efm_lengths);
    max_len = max(efm_lengths);
    
    title_str = sprintf('Distribution of EFM Lengths (N=%d)\nMean: %.1f | Median: %d | Range: [%d, %d]', ...
        size(E, 2), avg_len, med_len, min_len, max_len);
    title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
    
    % Add mean line
    xline(avg_len, '--r', 'LineWidth', 1.5, 'Label', 'Mean');
    
    % Save figure
    if nargin >= 2 && ~isempty(saveDir)
        if nargin >= 3
            filename = fullfile(saveDir, sprintf('%s_EFM_length_distribution.png', clusterName));
        else
            filename = fullfile(saveDir, 'EFM_length_distribution.png');
        end
        
        % Save as PNG and PDF for high quality
        saveas(fig, filename);
        saveas(fig, strrep(filename, '.png', '.pdf'));
        
        % Also save as .fig for MATLAB editing
        saveas(fig, strrep(filename, '.png', '.fig'));
        
       
    end
    
    close(fig);
end