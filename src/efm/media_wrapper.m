%% run_all_medias.m

% 0) Gurobi threads & paths (Do this once here)
nThreads = feature('numcores')-2;
fprintf('Using %d threads for Gurobi\n', nThreads);
setenv('OMP_NUM_THREADS', num2str(nThreads));
addpath(fullfile(getenv('GUROBI_HOME'),'examples','matlab'),'-begin');

% 1) Load irreversible pruned model ONCE to save I/O time
model_path = '/Users/benedikthaarscheidt/M.Sc./master_thesis_second_moment/Models/generic_models/E_coli/iML1515_pruned_permissive_biomass_loopless_v6.mat';
data = load(model_path, 'pruned_ir');
model_ir = data.pruned_ir;

% 2) Define Base Minimal Media Components
% These are the essential salts, minerals, and oxygen for aerobic E. coli growth
base_minimal = {'ex_nh4_e_b', 'ex_pi_e_b', 'ex_so4_e_b', 'ex_o2_e_b', ...
                'ex_h2o_e_b', 'ex_h_e_b', 'ex_k_e_b', 'ex_na1_e_b', ...
                'ex_mg2_e_b', 'ex_ca2_e_b', 'ex_cl_e_b', 'ex_fe2_e_b', ...
                'ex_fe3_e_b', 'ex_zn2_e_b', 'ex_mn2_e_b', 'ex_cu2_e_b', ...
                'ex_cobalt2_e_b', 'ex_mobd_e_b', 'ex_ni2_e_b', 'ex_sel_e_b', 'ex_tungs_e_b'};

% 3) Define Media Environments (Name, Active Carbon/Energy Sources)
% Expanded with common E. coli substrates using iML1515 terminology
media_configs = {
    'acetate',     [{'ex_ac_e_b'}, base_minimal];
    'fructose',    [{'ex_fru_e_b'}, base_minimal];
    'succinate',   [{'ex_succ_e_b'}, base_minimal];
    'pyruvate',    [{'ex_pyr_e_b'}, base_minimal];
    'lactate',     [{'ex_lac__D_e_b'}, base_minimal];
    'malate',      [{'ex_mal__L_e_b'}, base_minimal];
    'fumarate',    [{'ex_fum_e_b'}, base_minimal];
    'xylose',      [{'ex_xyl__D_e_b'}, base_minimal];
    'arabinose',   [{'ex_arab__L_e_b'}, base_minimal];
    'galactose',   [{'ex_gal_e_b'}, base_minimal];
};

% 4) Loop through and call the enumeration function
num_medias = size(media_configs, 1);
for i = 1:num_medias
    media_name = media_configs{i, 1};
    substrates = media_configs{i, 2};
    
    fprintf('\n======================================================\n');
    fprintf('STARTING WRAPPER FOR MEDIA: %s (%d/%d)\n', upper(media_name), i, num_medias);
    fprintf('======================================================\n');
    
    % Call the enumeration function with error handling
    try
        enumerate_efms_for_media(model_ir, media_name, substrates, nThreads);
    catch ME
        fprintf('\n[!] ERROR running %s: %s\n', upper(media_name), ME.message);
        fprintf('Skipping to the next media condition...\n');
        % It will automatically proceed to the next iteration
    end
end

fprintf('\n======================================================\n');
fprintf('All %d media enumerations completed!\n', num_medias);
fprintf('======================================================\n');