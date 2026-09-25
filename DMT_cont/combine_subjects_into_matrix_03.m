%% Step 3 of 4: Extracting variables from individual files (subject files) and saving them in a single matlab file
% for DMT_Project 1 data, both data sets: DMT_4cond and DMTcont

% Paths and parameters
% Input paths
%input_path ='...';
input_path = '...'; 

% Output paths
output_file = '/.../AllSubs_KOP_Curl_DMTcont_interpv2mil_LH.mat';
% output_file = '/.../AllSubs_KOP_Curl_DMTcont_interpv2mil_RH.mat';

%LH
conditions = {'DMT_all_left', 'PCB_all_left'}; 
% RH
% conditions = {'DMT_all_right', 'PCB_all_right'};

nConditions = 2;
nSubjects = 14;

%Preallocate cell arrays
all_KuraVertexGrid = cell(nConditions, nSubjects);
all_curlz = cell(nConditions, nSubjects);

% Loop through conditions and subjects
for cond = 1:nConditions
    for subj = 1:nSubjects
        fname = sprintf('Subj%02d_%s_RH.mat', subj, conditions{cond}); %Subj%02d_%s_LH.mat
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
save(output_file,'all_KuraVertexGrid', 'all_curlz', '-v7.3');

fprintf('Summary file saved: %s\n', output_file);
% Updated: 22Sep2026