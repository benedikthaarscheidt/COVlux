%% COVLUX MASTER PIPELINE ORCHESTRATOR
% =========================================================================
% Sequence: 
% 1. MRAS Mapping -> 2. Optimization -> 3. Master Comparison -> 
% 4. Biomass Distance -> 5. Biological Analysis -> 6. FBA Comparison
% =========================================================================


% Start a flight recorder to capture all command window output
diary('Pipeline_Run_Log.txt'); 

% 1. Define the relative paths to your scripts
pipeline_tasks = {
    %fullfile('src', 'mras', 'MRAS_mapping.m');
    fullfile('src', 'optimization', 'COVlux.m');
    fullfile('tests', 'master_comparison_cov.m');
    fullfile('tests', 'biomass_distance.m');
    fullfile('tests', 'analyze_biological_state.m');
    fullfile('tests', 'FBA_comp.m')
};

total_tasks = length(pipeline_tasks);
fprintf('Starting Pipeline at %s\n', datestr(now));
fprintf('-----------------------------------------------------------\n');

for i = 1:total_tasks
    [~, script_name, ~] = fileparts(pipeline_tasks{i});
    script_path = pipeline_tasks{i};
    
    fprintf('\n[%d/%d] RUNNING: %s\n', i, total_tasks, script_name);
    
    tic; 
    try
        % Using 'run' is essential here because the scripts are in 
        % different subdirectories. 'run' handles the context switching.
        run(script_path); 
        
        elapsed = toc;
        fprintf('[SUCCESS] Finished %s in %.2f seconds.\n', script_name, elapsed);
    catch ME
        elapsed = toc;
        fprintf('![CRASH] %s failed after %.2f seconds.\n', script_name, elapsed);
        fprintf('ERROR: %s\n', ME.message);
        
        % Log failure to a persistent file
        fid = fopen('pipeline_failures.txt', 'a');
        fprintf(fid, '[%s] Script: %s | Error: %s\n', datestr(now), script_name, ME.message);
        fclose(fid);
    end
    fprintf('-----------------------------------------------------------\n');
end

fprintf('\nPipeline Completed at %s\n', datestr(now));
diary off;