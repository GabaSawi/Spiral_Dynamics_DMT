%% Step 3 of 6: Computing Transition Probability Matrices for KL (Breaking of Detailed Balance)

%% Data and parameters
% clear all
% close all
addpath(genpath('/...')); %spm12
addpath(genpath('/.../1_Scripts/DMT_4cond/0_Obtaining_Raw_NMs/script_variabs&files')); %Glasser and parcels2schaefer
% Load flowmaps
addpath('/...'); % LH flowmaps
% addpath('/...'); % RH flowmaps

% load subjects mask
addpath('/.../0_Data/DMT_4cond')
load("mask_subjects_LH.mat")%LH
% load("mask_subjects_RH.mat") %RH

% check that mask loaded matches the hemisphere (example: condition 1, sub 1)
figure; imagesc(squeeze(mask_subs(1,1,:,:)))
figure; imagesc(squeeze(mask_subs2(1,1,:,:)))

downSRate = 2 ;                 % downsample the re-interpolation

% LH coordinates
xCord = -250:downSRate:250;
yCord = -150:downSRate:200;
posValid= gifti('Glasser360.L.flat.32k_fs_LR.surf.gii');

% RH coordinates
% xCord = -260:downSRate:240;
% yCord = -180:downSRate:170;
% posValid= gifti('Glasser360.R.flat.32k_fs_LR.surf.gii');

% Vertices, filter
NGROUP=4; % 1:pre-PCB / 2:post-PCB / 3:pre-DMT / 4:post-DMT

% Iterate over conditions
conditions = {'pre_PCB_left', 'post_PCB_left', 'pre_DMT_left', 'post_DMT_left'}; % uncomment for LH
% conditions = {'pre_PCB_right', 'post_PCB_right', 'pre_DMT_right', 'post_DMT_right'};  % uncomment for RH

NSUB=25;
NSUB_ok=17;
TR=2.00;  % Repetition Time (seconds)
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.01;                  % lowpass frequency of filter (Hz)
fhi = 0.08;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
%% Mask (create from valid vertices)
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
%% Map of neighbour indices kij2neigh(k,1:8) per grid-vertex-index kij(i,j)
Ni=size(yCord,2);   % Number rows
Nj=size(xCord,2);   % Number columns
Nneigh=8;           % Number neighbours per center-vertex
kij=reshape(1:(Ni*Nj),Nj,Ni)'; % CAREFUL transposed for L2R reading, NixNj changed

kij2neigh = zeros(Ni*Nj,Nneigh);
for i=1:Ni
    for j=1:Nj
        if mask(i,j)==1 % If center-bin within chicken
            k=kij(i,j);
            kij2neigh(k,1)=k+1;    % 1st (angle=0)
            kij2neigh(k,2)=k+Nj+1; % 2nd (angle=45)
            kij2neigh(k,3)=k+Nj;   % 3rd (angle=90)
            kij2neigh(k,4)=k+Nj-1; % 4th (angle=135)
            kij2neigh(k,5)=k-1;    % 5th (angle=180)
            kij2neigh(k,6)=k-Nj-1; % 6th (angle=225)
            kij2neigh(k,7)=k-Nj;   % 7th (angle=275)
            kij2neigh(k,8)=k-Nj+1; % 8th (angle=315)
        end
    end
end
fprintf('Per kij (vertex in grid list) neighbours index computed')

%% Loop over test subjects
for g=1:NGROUP 
    probtrans_group= zeros(Ni*Nj,Nneigh); 
    for sub= 1:NSUB_ok % WARNING: try and pass if not found (break for, continue)
        sub
        % VARIABLES TO SAVE (per subject!)
        probtrans_sub = NaN(Ni*Nj,Nneigh);
        % WARNING: if time removal, do it here (for now, keeping all tps)

        % Load Non-Normalized (Vx,Vy)
        output_path = '/...'; %lh
%         output_path = '/...'; %rh
        
        fname = sprintf('Subj%02d_%s_LH.mat', sub, conditions{g});
%         fname = sprintf('Subj%02d_%s_RH.mat', sub, conditions{g});
        load(fullfile(output_path, fname), 'vPhaseX', 'vPhaseY');
        % Check if file exists
        if exist("vPhaseX", 'var') == 1
            mask_SUB2=squeeze(mask_subs2(g,sub,:,:));
        else
            fprintf('File not found: %s. Skipping...\n', fname);
            continue;
        end

        Vx_flowmap_phase = vPhaseX;
        Vy_flowmap_phase = vPhaseY;
        clear vPhaseX vPhaseY;

        % Normalise the flowmaps
        Vx_flowmap_phase_norm= Vx_flowmap_phase./sqrt(Vx_flowmap_phase.^2 + Vy_flowmap_phase.^2);
        Vy_flowmap_phase_norm= Vy_flowmap_phase./sqrt(Vx_flowmap_phase.^2 + Vy_flowmap_phase.^2);

        % Variables needed: [0,2pi]-phase increments at center-vertex (Vx, Vy)->(moduleMatrix)
        phaseAngles = atan2(Vy_flowmap_phase_norm, Vx_flowmap_phase_norm); % Radians (-pi to pi)
        phaseAngles = mod(phaseAngles, 2*pi); % Wrap to [0, 2*pi]
        tps = size(Vx_flowmap_phase_norm,3); % Number of timepoints
        % Compute module at center-vertex
        moduleMatrix = sqrt(Vx_flowmap_phase_norm.^2 + Vy_flowmap_phase_norm.^2); % Magnitude (should be 1 if normalized)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Transition to neighbours at time t: 360º/8 BINS, ROTATED pi/4=45º divisions 
        probtrans = zeros(Ni*Nj,Nneigh); % (All vertices, neighbours)
        for t=1:tps %Sum all timepoints matrix to one matrix then divide by number of timepoints to average (instead of matrix per timepoint)
            moduleMatrix_t_mask=mask_SUB2.*moduleMatrix(:,:,t);
            for i=1:Ni
                for j=1:Nj 
                    if ~isnan(moduleMatrix_t_mask(i,j)) % Non-empty center vertex
                        % Counting neighbours counter-clockwise starting from angle=0
                        % 1st neighbour [2*pi-pi/8, 2*pi]U[0,pi/8]
                        if ((phaseAngles(i,j,t)>=(2*pi-pi/8) && phaseAngles(i,j,t)<=2*pi) | (phaseAngles(i,j,t)>=0 && (phaseAngles(i,j,t)<pi/8)))
                            if j~=Nj && ~isnan(moduleMatrix_t_mask(i,j+1))
                                probtrans(kij(i,j),1)=probtrans(kij(i,j),1)+moduleMatrix_t_mask(i,j);% 1st neighbour
                            end
                        else
                            % neighbours 2-to-8 --> Angles: [pi/8+(N-2)*pi/4, pi/8+(N-1)*pi/4] 2<=N<=8
                            for Neigh=2:8
                                if ((phaseAngles(i,j,t)>=(pi/8+(Neigh-2)*pi/4)) && (phaseAngles(i,j,t)<(pi/8+(Neigh-1)*pi/4)))
                                    if Neigh==2
                                        if i~=Ni && j~=Nj && ~isnan(moduleMatrix_t_mask(i+1,j+1))
                                            probtrans(kij(i,j),2)=probtrans(kij(i,j),2)+moduleMatrix_t_mask(i,j);% 2nd neighbour
                                        end
                                    elseif Neigh==3
                                        if i~=Ni && ~isnan(moduleMatrix_t_mask(i+1,j))
                                            probtrans(kij(i,j),3)=probtrans(kij(i,j),3)+moduleMatrix_t_mask(i,j);% 3rd neighhour
                                        end
                                    elseif Neigh==4
                                        if i~=Ni && j~=1 && ~isnan(moduleMatrix_t_mask(i+1,j-1))
                                            probtrans(kij(i,j),4)=probtrans(kij(i,j),4)+moduleMatrix_t_mask(i,j); % 4th neighbour
                                        end
                                    elseif Neigh==5
                                        if j~=1 && ~isnan(moduleMatrix_t_mask(i,j-1))
                                            probtrans(kij(i,j),5)=probtrans(kij(i,j),5)+moduleMatrix_t_mask(i,j); % 5th neighbour
                                        end
                                    elseif Neigh==6
                                        if i~=1 && j~=1 && ~isnan(moduleMatrix_t_mask(i-1,j-1))
                                            probtrans(kij(i,j),6)=probtrans(kij(i,j),6)+moduleMatrix_t_mask(i,j); % 6th neighbour
                                        end
                                    elseif Neigh==7
                                        if i~=1 && ~isnan(moduleMatrix_t_mask(i-1,j))
                                            probtrans(kij(i,j),7)=probtrans(kij(i,j),7)+moduleMatrix_t_mask(i,j); % 7th neighbour
                                        end
                                    elseif Neigh==8
                                        if i~=1 && j~=Nj && ~isnan(moduleMatrix_t_mask(i-1,j+1))
                                            probtrans(kij(i,j),8)=probtrans(kij(i,j),8)+moduleMatrix_t_mask(i,j); % 8th neighbour
                                        end
                                    end
                                end
                            end
                        end
                    end %non-NaN center
                end %for columns
            end %for rows
        end %for time points
        probtrans_sub=probtrans./tps; % Average

        output_path = '/.../DMT_4cond_ProbTrans_LH'; %lh
%         output_path = '/.../DMT_4cond_ProbTrans_RH'; %rh
        
        file_name_save = ['g', num2str(g), '_sub', num2str(sub), '_LH_probtrans.mat']; 
%         file_name_save = ['g', num2str(g), '_sub', num2str(sub), '_RH_probtrans.mat']; 
        full_path = fullfile(output_path, file_name_save);
        save(full_path, 'probtrans_sub');

        % GROUP AVERAGE PROBTRANS MATRIX
        probtrans_group = probtrans_group + probtrans_sub;
    end
    probtrans_group = probtrans_group./NSUB_ok;
    file_name_save=join(['G0',num2str(g),'_probtrans_LH.mat']); 
%     file_name_save=join(['G0',num2str(g),'_probtrans_RH.mat']);
    full_path = fullfile(output_path, file_name_save);
    save(full_path, 'probtrans_group');
end
% Updated: 22Sept2026