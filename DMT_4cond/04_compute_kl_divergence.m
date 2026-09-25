%% Step 4 of 6: Computing KL (Breaking of Detailed Balance)

%% Data and parameters
addpath(genpath('/...'));%spm12
addpath(genpath('/.../1_Scripts/DMT_4cond/0_Obtaining_Raw_NMs/script_variabs&files'));

% set input and output paths
% LH
addpath('/.../0_Data/DMT_4cond/DMT_4cond_LH_20Jan26/DMT_4cond_ProbTrans_LH')
output_path = '/.../0_Data/DMT_4cond/DMT_4cond_LH_20Jan26/DMT_4cond_KLfromProbTrans_LH';
% RH
% addpath('/.../0_Data/DMT_4cond/DMT_4cond_RH_20Jan26/DMT_4cond_ProbTrans_RH');
% output_path = '/.../0_Data/DMT_4cond/DMT_4cond_RH_20Jan26/DMT_4cond_KLfromProbTrans_RH';

NGROUP=4; % 1:pre-PCB / 2:post-PCB / 3:pre-DMT / 4:post-DMT
NSUB=25;
NSUB_ok=17;

%% Grid variables
downSRate = 2 ;                        % downsample the re-interpolation

% LH coordinates
xCord = -250:downSRate:250;
yCord = -150:downSRate:200;
posValid= gifti('Glasser360.L.flat.32k_fs_LR.surf.gii');
% RH coordinates
% xCord = -260:downSRate:240;
% yCord = -180:downSRate:170;
% posValid= gifti('Glasser360.R.flat.32k_fs_LR.surf.gii');

% Indices over grid kij(i,j)
N = size(posValid.vertices,1);
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

% k-index on grid vertices
Nneigh=8; % Number neighbours per center-vertex
kij=reshape(1:(Ni*Nj),Nj,Ni)'; % CAREFUL transposed for L2R reading, NixNj changed

% Interpolate Schaefer 1000 labels onto regular grid + masking 
DLfile = 'Schaefer2018_1000Parcels_7Networks_order.dscalar.nii';
fslr_7os= ft_read_cifti(DLfile);                                 % Schaefer labels (per grayordinate)
vertex2mil = fslr_7os.dscalar(1:N); % LH
% vertex2mil = fslr_7os.dscalar(N+1:end); % RH

% interpolate vertex2mil (from github 7Networks) % Sanitized interp_v2mil_Yeo7 (from github Yeo7)
interp_v2mil_Yeo7(:,:) = griddata(x,y,vertex2mil,X, Y,'nearest') ;
interp_v2mil_Yeo7 = interp_v2mil_Yeo7.*mask ;


% Variables involved in computing KL
ratio_prob_inf_subs = zeros(NGROUP,NSUB_ok,Ni*Nj,Nneigh); % k(neigh): P(neigh->k)=0
KL_prob_subs = nan(NGROUP,NSUB_ok,Ni*Nj);
KL_prob_grid_subs = nan(NGROUP,NSUB_ok,Ni,Nj);


%% 1. Compute probability of phase flow directions between neighbours (ratio forward/reverse flow)
eps_Laplace = 1e-10; %
for g=1:NGROUP
    nsub=0;
    for sub= 1:NSUB_ok % WARNING: try and pass if not found (break for, continue)
        sub
        file_name_open=['g', num2str(g), '_sub', num2str(sub), '_LH_probtrans.mat'];
%         file_name_open=['g', num2str(g), '_sub', num2str(sub), '_RH_probtrans.mat'];
        % Check if file exists
        if exist(file_name_open, 'file') == 2
            load(file_name_open)
            nsub=nsub+1;
        else
            fprintf('File not found: %s. Skipping...\n', file_name_open);
            continue;
        end

        probtrans_sub=probtrans_sub+eps_Laplace;
        
        % Ratio per pair of (vertex,neigh) Pa->b/Pb->a
        ratioprobtrans = nan(Ni*Nj,Nneigh);
        ratioprobtrans_inf = zeros(Ni*Nj,Nneigh);

        % 1.1.1. Get ratio Pab/Pba, for later Pab*log(Pab/Pba)
        for i=1:Ni
            for j=1:Nj
                if mask(i,j)==1 % Vertex in mask
                    k=kij(i,j);
                    if mask(i,j+1)==1 & probtrans_sub(k+1,5)~=0 % 1st Neighbour in mask
                        ratioprobtrans(k,1)=probtrans_sub(k,1)/probtrans_sub(k+1,5);
                        if probtrans_sub(k+1,5)<eps_Laplace; ratioprobtrans_inf(k,1)=1; end % Count infinite P(neigh->k)
                    end
                    if mask(i+1,j+1)==1 & probtrans_sub(k+Nj+1,6)~=0 % 2nd Neighbour in mask
                        ratioprobtrans(k,2)=probtrans_sub(k,2)/probtrans_sub(k+Nj+1,6);
                        if probtrans_sub(k+Nj+1,6)<eps_Laplace; ratioprobtrans_inf(k,2)=1; end
                    end
                    if mask(i+1,j)==1 & probtrans_sub(k+Nj,7)~=0 % 3rd Neighbour in mask
                        ratioprobtrans(k,3)=probtrans_sub(k,3)/probtrans_sub(k+Nj,7);
                        if probtrans_sub(k+Nj,7)<eps_Laplace; ratioprobtrans_inf(k,3)=1; end
                    end
                    if mask(i+1,j-1)==1 & probtrans_sub(k+Nj-1,8)~=0 % 4th Neighbour in mask
                        ratioprobtrans(k,4)=probtrans_sub(k,4)/probtrans_sub(k+Nj-1,8);
                        if probtrans_sub(k+Nj-1,8)<eps_Laplace; ratioprobtrans_inf(k,4)=1; end
                    end
                    if mask(i,j-1)==1 & probtrans_sub(k-1,1)~=0 % 5th Neighbour in mask
                        ratioprobtrans(k,5)=probtrans_sub(k,5)/probtrans_sub(k-1,1);
                        if probtrans_sub(k-1,1)<eps_Laplace; ratioprobtrans_inf(k,5)=1; end
                    end
                    if mask(i-1,j-1)==1 & probtrans_sub(k-Nj-1,2)~=0 % 6th Neighbour in mask
                        ratioprobtrans(k,6)=probtrans_sub(k,6)/probtrans_sub(k-Nj-1,2);
                        if probtrans_sub(k-Nj-1,2)<eps_Laplace; ratioprobtrans_inf(k,6)=1; end
                    end
                    if mask(i-1,j)==1 & probtrans_sub(k-Nj,3)~=0 % 7th Neighbour in mask
                        ratioprobtrans(k,7)=probtrans_sub(k,7)/probtrans_sub(k-Nj,3);
                        if probtrans_sub(k-Nj,3)<eps_Laplace; ratioprobtrans_inf(k,7)=1; end
                    end
                    if mask(i-1,j+1)==1 & probtrans_sub(k-Nj+1,4)~=0 % 8th Neighbour in mask
                        ratioprobtrans(k,8)=probtrans_sub(k,8)/probtrans_sub(k-Nj+1,4);
                        if probtrans_sub(k-Nj+1,4)<eps_Laplace; ratioprobtrans_inf(k,8)=1; end
                    end
                end
            end
        end
        ratio_prob_inf_subs(g,nsub,:,:)=ratioprobtrans_inf; % save count of infinite ratio prob trans
        
        % COMPUTE KL per vertex
        for k=1:Ni*Nj
            PijPji = ratioprobtrans(k,:); 
            kl_prob(k)=sum(probtrans_sub(k,:).*log(PijPji));
        end
        KL_prob_subs(g,nsub,:)=kl_prob; % SAVE KL per vertex
        disp(['KL g:',num2str(g),' sub:',num2str(sub),' computed'])
        
        file_name_save=join(['g',num2str(g),'sub', num2str(sub),'LH_KL.mat']);
        %file_name_save=join(['g',num2str(g),'sub', num2str(sub),'RH_KL.mat']);
        full_path = fullfile(output_path,file_name_save);
        save(full_path,'kl_prob');
        % KL per vertex -> GRID
        kl_prob_grid = reshape(kl_prob,Nj,Ni)'; % Reshape to GRID
        KL_prob_grid_subs(g,nsub,:,:)=kl_prob_grid;
    end
end
file_name_save=join(['All_groups_subs_LH','_KL.mat']);
% file_name_save=join(['All_groups_subs_RH','_KL.mat']);
full_path = fullfile(output_path,file_name_save);
%save(full_path,'KL_prob_grid_subs','KL_prob_subs','KL_prob_grid_subs');

% Updated: 22Sept2026