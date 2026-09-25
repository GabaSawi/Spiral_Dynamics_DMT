%% Step 6 of 6: Computing MEAN and STD of KOP and Curl for Kuramoto and Spiral Vorticity and Turbulence
% Pipeline order: 01_compute_kop_curl_phaseflow -> 02_build_subject_masks -> 03_compute_transition_probabilities
%              -> 04_compute_kl_divergence -> 05_combine_subjects_into_matrix -> 06 (this file)
% Output of this script (DMT_4cond_allNMs_newYeo7_RH_vf.mat) is what the figure scripts load directly.
% Define NMs:
NM = {'KOP','Curl','KOP SD','Curl SD','KL'};

% Add paths:
addpath(genpath('/...')); % SPM12 path
addpath(genpath('/.../1_Scripts/DMT_4cond/0_Obtaining_Raw_NMs/script_variabs&files'));

% Add paths to data
% KOP and Curl
addpath('/.../0_Data/DMT_4cond/DMT_4cond_LH_20Jan26');%LH
% addpath('/.../0_Data/DMT_4cond/DMT_4cond_RH_20Jan26'); %RH

% KL
addpath('.../0_Data/DMT_4cond/DMT_4cond_LH_20Jan26/DMT_4cond_KLfromProbTrans_LH'); %LH
% addpath('/.../0_Data/DMT_4cond/DMT_4cond_RH_20Jan26/DMT_4cond_KLfromProbTrans_RH') %RH
% % % addpath('/home/user-cbc/Desktop/PhD/Projects/Turbulence_Spirals/Results/Results_KL_final') % WARNING: DMTcont KL 

% loading data LH data
load("AllSubs_KOP_Curl_DMT_4cond_interpv2mil_LH.mat", 'all_KuraVertexGrid', 'all_curlz'); % Grid results were not affected by old grid2parcels nor vertex2mil
load('All_groups_subs_LH_KL.mat','KL_prob_grid_subs'); 
% RH data
% load("AllSubs_KOP_Curl_DMT_4cond_interpv2mil_RH.mat", 'all_KuraVertexGrid', 'all_curlz'); % Grid results were not affected by old grid2parcels nor vertex2mil
% load('All_groups_subs_RH_KL.mat','KL_prob_grid_subs'); 

% % % Check 4cond 17subj
% for cond=1:4
%     figure
%     for subj=1:17
%         subplot(3,6,subj)
%         imagesc(squeeze(KL_prob_grid_subs(cond,subj,:,:)))
%     end
% end
 
% Preallocate arrays
all_Kura1000_mean = cell(4,17);
all_spiral1000_mean = cell(4,17);
all_KL1000_mean = cell(4,17);
all_KuraGrid1000_STD = cell(4,17);
all_spiral1000_STD = cell(4,17);

%% Processing the data
downSRate = 2;                        % Downsample the re-interpolation
% LH coordinates
xCord = -250:downSRate:250;
yCord = -150:downSRate:200;
posValid = gifti('Glasser360.L.flat.32k_fs_LR.surf.gii');
% RH coordinates
% xCord = -260:downSRate:240; 
% yCord = -180:downSRate:170;
% posValid = gifti('Glasser360.R.flat.32k_fs_LR.surf.gii');

N = size(posValid.vertices,1);
x = double(posValid.vertices(:, 1));
y = double(posValid.vertices(:, 2));
k_shape = alphaShape(x, y, 4);
[a, b] = k_shape.boundaryFacets();
bw = poly2mask(b(:,1)-min(xCord)+1, b(:,2)-min(yCord)+1, max(yCord)-min(yCord)+1, max(xCord)-min(xCord)+1);
mask = double(bw(1:downSRate:end, 1:downSRate:end));
mask(mask == 0) = nan;
[X_grid, Y_grid] = meshgrid(xCord, yCord);

% Instead of grid2parcels use (Yeo7 from github) vertex2mil --> interpolate + masking 
addpath('/media/user-cbc/Data/PhD/Projects/DMT_Experiece_Tracking_4GitHub/0_Data')
DLfile = 'Schaefer2018_1000Parcels_7Networks_order.dscalar.nii';
fslr_7os= ft_read_cifti(DLfile); % Schaefer labels (per grayordinate)
% LH
%vertex2mil = fslr_7os.dscalar(1:N);
%RH
vertex2mil = fslr_7os.dscalar(N+1:end);

% interpolate vertex2mil (from github 7Networks) % Sanitized interp_v2mil_Yeo7 (from github Yeo7)
interp_v2mil_Yeo7(:,:) = griddata(x,y,vertex2mil,X_grid, Y_grid,'nearest') ;
interp_v2mil_Yeo7 = interp_v2mil_Yeo7.*mask ;

%% Iterate over data
% LH
seedsLH = 1:500;
% RH 
% seedsRH = 501:1000;

nSeed = numel(seedsRH);          % 500
for cond = 1:size(all_KuraVertexGrid, 1)
    for subj = 1:size(all_KuraVertexGrid, 2)
        
        fprintf('  Subject %d\n', subj);
        
        grid_kop_sub = all_KuraVertexGrid{cond,subj};
        grid_curl_sub = all_curlz{cond,subj};
        grid_kl_sub = squeeze(KL_prob_grid_subs(cond,subj,:,:));
        
        T = size(grid_kop_sub,3);

        % Subject variables to compute
        KuraGrid1000_mean = nan(nSeed, T);
        spiral1000_mean   = nan(nSeed, T);
        KuraGrid1000_STD  = nan(nSeed, T);
        spiral1000_STD    = nan(nSeed, T);
        KL1000_mean       = nan(nSeed, 1);

        % From NM in vertices --> 1000 (500) Parcellation  
        % Kuramoto & Curl from interp_v2mil_Yeo7
        for si = 1:nSeed
            seed = seedsRH(si);
            mask_roi = (interp_v2mil_Yeo7 == seed);                 % logical [Ny x Nx]
            % reshape to pixels x time
            K = reshape(grid_kop_sub, [], T);                       % [Ny*Nx x T]
            C = reshape(grid_curl_sub, [], T);            
            KL= grid_kl_sub(:);
        
            idx = find(mask_roi(:));                           % linear indices of this seed’s pixels
       
            % Mean across pixels (omit NaNs)
            KuraGrid1000_mean(si,:) = mean(abs(K(idx,:)),1,'omitnan');
            spiral1000_mean(si,:)   = mean(abs(C(idx,:)),1,'omitnan');
            KL1000_mean(si)   = mean(abs(KL(idx,:)),1,'omitnan');
            
            % STD across pixels (omit NaNs)
            KuraGrid1000_STD(si,:) = std((K(idx,:)),0, 1,'omitnan');
            spiral1000_STD(si,:) = std(abs(C(idx,:)),0, 1,'omitnan');
        end

        all_Kura1000_mean{cond,subj} = KuraGrid1000_mean; 
        all_spiral1000_mean{cond,subj} = spiral1000_mean;
        all_KL1000_mean{cond,subj} = KL1000_mean;
        all_KuraGrid1000_STD{cond,subj} = KuraGrid1000_STD;
        all_spiral1000_STD{cond,subj} = spiral1000_STD;
 
        clear KuraGrid1000_STD spiral1000_STD KuraGrid1000_mean spiral1000_mean KL1000_mean K C KL;
       
    end
end

%% Save
%LH
save('/.../0_Data/DMT_4cond/DMT_4cond_LH_20Jan26/DMT_4cond_allNMs_newYeo7_LH_vf.mat','all_Kura1000_mean','all_spiral1000_mean','all_KL1000_mean', 'all_KuraGrid1000_STD', 'all_spiral1000_STD');
%RH
% save('/.../0_Data/DMT_4cond/DMT_4cond_RH_20Jan26/DMT_4cond_allNMs_newYeo7_RH_vf.mat','all_Kura1000_mean','all_spiral1000_mean','all_KL1000_mean', 'all_KuraGrid1000_STD', 'all_spiral1000_STD');
fprintf('Data Saved Successfully!')
% Updated: 20Feb2026