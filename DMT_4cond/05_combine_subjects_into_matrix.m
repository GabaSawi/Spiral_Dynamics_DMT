%% Step 5 of 6: Extracting variables from individual files (subject files) and saving them in a single matlab file
% Pipeline order: 01_compute_kop_curl_phaseflow -> 02_build_subject_masks -> 03_compute_transition_probabilities
%              -> 04_compute_kl_divergence -> 05 (this file) -> 06_aggregate_yeo7_parcels
% Paths and parameters

% Paths for input data
%input_path = '/media/user-cbc/Data/PhD/Projects/DMT_Experiece_Tracking_4GitHub/0_Data/DMT_4cond/DMT_4cond_LH_20Jan26'; %LH
input_path = '/media/user-cbc/Data/PhD/Projects/DMT_Experiece_Tracking_4GitHub/0_Data/DMT_4cond/DMT_4cond_RH_20Jan26'; %rh

% Output paths
%output_file ='/media/user-cbc/Data/PhD/Projects/DMT_Experiece_Tracking_4GitHub/0_Data/DMT_4cond/DMT_4cond_LH_20Jan26/AllSubs_KOP_Curl_DMT_4cond_interpv2mil_LH.mat'; %LH
output_file = '/media/user-cbc/Data/PhD/Projects/DMT_Experiece_Tracking_4GitHub/0_Data/DMT_4cond/DMT_4cond_RH_20Jan26/AllSubs_KOP_Curl_DMT_4cond_interpv2mil_RH.mat'; %RH

% Define conditions
%conditions = {'pre_PCB_left', 'post_PCB_left', 'pre_DMT_left','post_DMT_left'}; %LH
conditions = {'pre_PCB_right', 'post_PCB_right', 'pre_DMT_right', 'post_DMT_right'}; %RH

nConditions = 4;
nSubjects = 17;
%% Preallocate cell arrays
all_KuraVertexGrid = cell(nConditions, nSubjects);
all_curlz = cell(nConditions, nSubjects);

% Loop through conditions and subjects
for cond = 1:nConditions
    for subj = 1:nSubjects
        %fname = sprintf('Subj%02d_%s_LH.mat', subj, conditions{cond}); %lh
        fname = sprintf('Subj%02d_%s_RH.mat', subj, conditions{cond}); %rh
        
        fpath = fullfile(input_path, fname);

        if isfile(fpath)
            S = load(fpath,'KuraVertexGrid','curlz');
            all_KuraVertexGrid{cond, subj}     = S.KuraVertexGrid;
            all_curlz{cond, subj}              = S.curlz;

        else
            warning('Missing file: %s', fpath);
        end
    end
end

% Save to single file
%save(output_file,'all_KuraVertexGrid', 'all_curlz', '-v7.3');

fprintf('Summary file saved: %s\n', output_file);
%Updated: 20Feb2026
