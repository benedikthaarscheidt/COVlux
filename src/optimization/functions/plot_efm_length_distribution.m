function plot_efm_length_distribution(E, rxnNames, saveDir, clusterName)
    
    % Calculate lengths
    efm_lengths = sum(abs(E) > 1e-9, 1);
    
    % Create Plot (Made taller to fit the table)
    fig = figure('Name', 'EFM Length Analysis', 'Color', 'w', 'Visible', 'off', 'Position', [100, 100, 900, 800]);
    
    % --- TOP HALF: HISTOGRAM ---
    subplot(2, 1, 1);
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
    
    % --- BOTTOM HALF: TABLE OF SHORTEST EFMs ---
    % Find the 10 shortest EFMs
    [sorted_lengths, sort_idx] = sort(efm_lengths, 'ascend');
    num_to_show = min(10, length(efm_lengths));
    shortest_idx = sort_idx(1:num_to_show);
    
    % Prepare table data
    table_data = cell(num_to_show, 3);
    for i = 1:num_to_show
        e_idx = shortest_idx(i);
        len = sorted_lengths(i);
        
        % Extract reaction names for this EFM
        active_rxns = rxnNames(abs(E(:, e_idx)) > 1e-9);
        rxn_str = strjoin(active_rxns, ', ');
        
        % Truncate string if it's too long so it fits on the PDF page
        if length(rxn_str) > 140
            rxn_str = [rxn_str(1:137), '...'];
        end
        
        table_data{i, 1} = sprintf('EFM %d', e_idx);
        table_data{i, 2} = len;
        table_data{i, 3} = rxn_str;
    end
    
    % Create table UI component
    ax2 = subplot(2, 1, 2);
    axis(ax2, 'off'); % Hide the axes, we just want the table
    
    uitable(fig, 'Data', table_data, ...
        'ColumnName', {'Original Index', 'Length', 'Active Reactions'}, ...
        'Units', 'normalized', ...
        'Position', [0.05 0.05 0.9 0.4], ... % Position in the bottom half
        'ColumnWidth', {100, 60, 680}, ...   % Width constraints
        'RowName', [], ...
        'FontSize', 9);
    
    % --- SAVE OUTPUTS ---
    if nargin >= 3 && ~isempty(saveDir)
        if nargin >= 4
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