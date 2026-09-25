%% Step 4 of 4: Computing MEAN and STD of KOP and Curl for Kuramoto and Spiral Vorticity and Turbulence
% Output of this script (DMTcont_allNMs_newYeo7_LH_vf.mat) is what the figure scripts load directly.

% Define NMs:
NM = {'KOP','Curl','KOP SD','Curl SD'};

% Add paths:
addpath(genpath('/...')); % SPM12 path
addpath(genpath('/.../1_Scripts/DMT_cont/0_Obtaining_Raw_NMs/script_variabs&files')); %functions path

% loading data for KOP, Curl
%LH
load("AllSubs_KOP_Curl_DMTcont_interpv2mil_LH.mat", 'all_KuraVertexGrid', 'all_curlz');

% RH
% load("AllSubs_KOP_Curl_DMTcont_interpv2mil_RH.mat", 'all_KuraVertexGrid', 'all_curlz');

% Saving variables
all_Kura1000_mean = cell(2,14);
all_spiral1000_mean = cell(2,14);
all_KuraGrid1000_STD = cell(2,14);
all_spiral1000_STD = cell(2,14);

%% Processing the data
downSRate = 2; 
% LH coordinates
xCord = -250:downSRate:250;
yCord = -150:downSRate:200;
posValid= gifti('Glasser360.L.flat.32k_fs_LR.surf.gii');

% RH coordinates
% xCord = -260:downSRate:240;
% yCord = -180:downSRate:170;
% posValid= gifti('Glasser360.R.flat.32k_fs_LR.surf.gii');

% Indices over grid kij(i,j)
Ni=size(yCord,2);   % Number rows
Nj=size(xCord,2);   % Number columns
x = double(posValid.vertices(:,1));
y = double(posValid.vertices(:,2));

% k = convhull(x,y);
k = alphaShape(x,y,4) ;
[a, b] = k.boundaryFacets();
bw = poly2mask(b(:,1)-min(xCord)+1,b(:,2)-min(yCord)+1,max(yCord)-min(yCord)+1,...
    max(xCord)-min(xCord)+1);
mask = double(bw(1:downSRate:end,1:downSRate:end)) ;
mask(mask==0) = nan ;
[X,Y] = meshgrid(xCord, yCord);
fprintf('Mask created')

N = size(posValid.vertices,1);

% Interpolate Schaefer 1000 labels onto regular grid + masking 
DLfile = 'Schaefer2018_1000Parcels_7Networks_order.dscalar.nii';
fslr_7os= ft_read_cifti(DLfile);                                 % Schaefer labels (per grayordinate)

vertex2mil = fslr_7os.dscalar(1:N); % LH
% vertex2mil = fslr_7os.dscalar(N+1:end); % RH

% interpolate vertex2mil (from github 7Networks) % Sanitized interp_v2mil_Yeo7 (from github Yeo7)
interp_v2mil_Yeo7(:,:) = griddata(x,y,vertex2mil,X, Y,'nearest') ;
interp_v2mil_Yeo7 = interp_v2mil_Yeo7.*mask ;

%% Iterate over data
% LH
%seedsLH = 1:500;
% RH 
seedsRH = 501:1000;

nSeed = numel(seedsRH);          % 500
for cond = 1:size(all_KuraVertexGrid, 1)
    for subj = 1:size(all_KuraVertexGrid, 2)
        
        fprintf('  Subject %d\n', subj);
        
        grid_kop_sub = all_KuraVertexGrid{cond,subj};
        grid_curl_sub = all_curlz{cond,subj};

        T = size(grid_kop_sub,3);

        % Subject variables to compute
        KuraGrid1000_mean = nan(nSeed, T);
        spiral1000_mean   = nan(nSeed, T);
        KuraGrid1000_STD  = nan(nSeed, T);
        spiral1000_STD    = nan(nSeed, T);

        % From NM in vertices --> 1000 (500) Parcellation
        % Kuramoto & Curl from interp_v2mil_Yeo7
        for si = 1:nSeed
            seed = seedsRH(si);
            mask_roi = (interp_v2mil_Yeo7 == seed);                 % logical [Ny x Nx]
            % reshape to pixels x time
            K = reshape(grid_kop_sub, [], T);                       % [Ny*Nx x T]
            C = reshape(grid_curl_sub, [], T);

            idx = find(mask_roi(:));                           % linear indices of this seed’s pixels

            % Mean across pixels (omit NaNs)
            KuraGrid1000_mean(si,:) = mean(abs(K(idx,:)),1,'omitnan');
            spiral1000_mean(si,:)   = mean(abs(C(idx,:)),1,'omitnan');

            % STD across pixels (omit NaNs)
            KuraGrid1000_STD(si,:) = std((K(idx,:)),0, 1,'omitnan');
            spiral1000_STD(si,:) = std(abs(C(idx,:)),0, 1,'omitnan');
        end

        % Store
        all_Kura1000_mean{cond,subj} = KuraGrid1000_mean;
        all_spiral1000_mean{cond,subj} = spiral1000_mean;
        all_KuraGrid1000_STD{cond,subj} = KuraGrid1000_STD;
        all_spiral1000_STD{cond,subj} = spiral1000_STD;

        clear KuraGrid1000_STD spiral1000_STD KuraGrid1000_mean spiral1000_mean K C;
       
    end
end

%% Save
% LH
save('/.../0_Data/DMT_cont/DMT_cont_LH_20Jan26/DMTcont_allNMs_newYeo7_LH_vf.mat','all_Kura1000_mean','all_spiral1000_mean', 'all_KuraGrid1000_STD', 'all_spiral1000_STD');
%RH
% save('/.../0_Data/DMT_cont/DMT_cont_RH_20Jan26/DMTcont_allNMs_newYeo7_RH_vf.mat','all_Kura1000_mean','all_spiral1000_mean', 'all_KuraGrid1000_STD', 'all_spiral1000_STD');
fprintf('Data Saved Successfully!')
% Updated: 22Sep2026