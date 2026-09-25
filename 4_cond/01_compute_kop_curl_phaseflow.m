%% Step 1 of 6: Computing KOP and Curl for Kuramoto and Spiral Vorticity and Turbulence

clear all
close all

% Adding the paths.
addpath(genpath('/...')); % SPM12 path
addpath(genpath('/.../1_Scripts/DMT_4cond/0_Obtaining_Raw_NMs/script_variabs&files')); %Glasser and parcels2schaefer
addpath(genpath('/.../DMT_HCP_32k')); % Raw data matrix

% Set output path
output_path = '/...'; % uncomment for LH

% Loading variables and setting coordinates for the grid
downSRate = 2;                        % Downsample the re-interpolation

% LH coordinates
xCord = -250:downSRate:250;
yCord = -150:downSRate:200;
posValid = gifti('Glasser360.L.flat.32k_fs_LR.surf.gii');
posValid_xy=double(posValid.vertices(:, 1:2));
% RH coordinates
% xCord = -260:downSRate:240; 
% yCord = -180:downSRate:170;
% posValid = gifti('Glasser360.R.flat.32k_fs_LR.surf.gii');
posValid_xy=double(posValid.vertices(:, 1:2));

% Loading project data
% LH data
load("DMT_project_data_left.mat");
% RH data
% load('DMT_project_data_right.mat');

% Setting the frequency filter
TR = 2.00;                            % Repetition Time (seconds)
fnq = 1/(2*TR);                       % Nyquist frequency
flp = 0.01;                          % Lowpass frequency of filter (Hz)
fhi = 0.08;                           % Highpass frequency
Wn = [flp/fnq fhi/fnq];               % Butterworth bandpass non-dimensional frequency
k = 2;                                % 2nd order Butterworth filter
[bfilt, afilt] = butter(k, Wn);       % Construct the filter

% Iterate over conditions
% LH
conditions = {'pre_PCB_left', 'post_PCB_left', 'pre_DMT_left','post_DMT_left'}; %RH: comment
% RH
% conditions = {'pre_PCB_right', 'post_PCB_right', 'pre_DMT_right', 'post_DMT_right'};

for cond = 1:length(conditions)
    fprintf('Processing Condition: %s\n', conditions{cond});
    data = evalin('base', conditions{cond});

    for subj = 1:length(data)
        fprintf('  Subject %d\n', subj);
        ts = data{subj, 1};
        N = size(ts, 1);
        tss = ts;
        T = size(tss,2);

        % Temporal Smoothing (replacing zeros with NaNs for global measures)
        signal_filt = zeros(N, T);
        for seed = 1:N
            tss(seed,:) = detrend(tss(seed,:) - nanmean(tss(seed,:)));
            if ts(seed,1) == 0
                signal_filt(seed,:) = nan;
            else
                signal_filt(seed, :) = filtfilt(bfilt, afilt, tss(seed, :)); %temporal smoothing
            end
        end

        % Spatial Smoothing
        dist = zeros(1, N);
        ts_filtered = zeros(N, T);
        for i = 1:N
            membervec=[];
            for j = 1:N
                dist(j) = sum((posValid_xy(i, :) - posValid_xy(j, :)).^2);
            end
            [aux, index] = sort(dist, 'ascend');
            membervec = index(1:60);
            gaussweight = exp(-aux(1:60)/16);
            ts_filtered(i, :) = gaussweight * signal_filt(membervec, :);
        end

        % Compute phases for Kuramoto
        Phases = zeros(N, T);
        for seed = 1:N
            ts_filtered(seed, :) = detrend(ts_filtered(seed, :) - nanmean(ts_filtered(seed, :)));
            Xanalytic = hilbert(ts_filtered(seed, :)); % mean subtraction redundant previously
            Phases(seed, :) = angle(Xanalytic);
        end

        % Create a Regular Mask
        dataIn = ts_filtered;
        x = double(posValid.vertices(:, 1));
        y = double(posValid.vertices(:, 2));
        k_shape = alphaShape(x, y, 4);
        [a, b] = k_shape.boundaryFacets();

        bw = poly2mask(b(:,1)-min(xCord)+1, b(:,2)-min(yCord)+1, max(yCord)-min(yCord)+1, max(xCord)-min(xCord)+1);

        mask = double(bw(1:downSRate:end, 1:downSRate:end));
        mask(mask == 0) = nan;

        [X_grid, Y_grid] = meshgrid(xCord, yCord);
        
        % Interpolate the data
        sigSSpec = dataIn;
        sigSSpec_reshape = nan(length(yCord), length(xCord), size(dataIn, 2));
        for iTime = 1:size(dataIn, 2)
            sigSSpec_reshape(:, :, iTime) = griddata(x, y, sigSSpec(:, iTime), X_grid, Y_grid, 'cubic');
            sigSSpec_reshape(:, :, iTime) = sigSSpec_reshape(:, :, iTime) .* mask;
        end
        clear dataIn

        %% Compute the Phase Vector Field, KOP, and Curl
        temp1_smooth_phase_sur = nan(size(sigSSpec_reshape));
        for irow = 1:size(temp1_smooth_phase_sur, 1)
            for icol = 1:size(temp1_smooth_phase_sur, 2)
                temp1 = sigSSpec_reshape(irow,icol,:);
                if nansum(abs(temp1(:))) ~= 0
                    temp1_smooth_phase_sur(irow,icol,:) = angle(hilbert(temp1(:)));
                end
            end
        end

        phaseSig = temp1_smooth_phase_sur;
        Vx_flowmap_norm_phase = [];
        Vy_flowmap_norm_phase = [];
        vPhaseX = zeros(size(phaseSig));
        vPhaseY = zeros(size(phaseSig));
        KuraVertexGrid = nan(size(phaseSig));

        for iTime = 1:size(phaseSig,3)
            for iX = 1:size(phaseSig,1)
                vPhaseX(iX,2:end-1,iTime) = (anglesubtract(phaseSig(iX,3:end,iTime), phaseSig(iX,1:end-2,iTime)))/2;
            end
            for iY = 1:size(phaseSig,2)
                vPhaseY(2:end-1,iY,iTime) = (anglesubtract(phaseSig(3:end,iY,iTime), phaseSig(1:end-2,iY,iTime)))/2;
            end
            for iX = 2:size(phaseSig,1)-1
                for iY = 2:size(phaseSig,2)-1
                    effNum = 9 - length(find(isnan(phaseSig(iX-1:iX+1,iY-1:iY+1,iTime))));
                    KuraVertexGrid(iX,iY,iTime) = abs(nansum(nansum(exp(1i*phaseSig(iX-1:iX+1,iY-1:iY+1,iTime)))))/effNum;
                end
            end
        end

        Vx_flowmap_norm_phase = -vPhaseX ./ sqrt(vPhaseX.^2 + vPhaseY.^2);
        Vy_flowmap_norm_phase = -vPhaseY ./ sqrt(vPhaseX.^2 + vPhaseY.^2);
        %clear vPhaseX vPhaseY

        % Compute Curl
        curlz = [];
        cav = [];
        for time = 1:size(Vx_flowmap_norm_phase, 3)
            temp1_vx = Vx_flowmap_norm_phase(:,:,time);
            temp1_vy = Vy_flowmap_norm_phase(:,:,time);
            [x_curl, y_curl] = meshgrid(1:size(temp1_vx, 2), 1:size(temp1_vx, 1));
            [curlz(:,:,time), cav(:,:,time)] = curl(x_curl, y_curl, temp1_vx, temp1_vy);
        end


        % Save subject results
        % LH
        fname = sprintf('Subj%02d_%s_LH.mat', subj, conditions{cond}); 
        % RH
%         fname = sprintf('Subj%02d_%s_RH.mat', subj, conditions{cond}); 
        save(fullfile(output_path, fname), 'KuraVertexGrid', 'curlz', 'vPhaseX', 'vPhaseY','-v7.3');

        clearvars KuraVertexGrid curlz cav temp1_smooth_phase_sur sigSSpec_reshape Vx_flowmap_norm_phase Vy_flowmap_norm_phase Phases ts ts_filtered tss

    end
end
% Updated: 22Sept2026