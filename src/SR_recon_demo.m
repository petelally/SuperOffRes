% SR_recon_demo takes the provided Bruker data and performs multi-freq 
% super-res reconstruction as outlined in [REF]. 
%
% For each frequency offset, a super-res image is reconstructed:
% 1) Constructs matrix of k-space bands and pattern shifts (psi)
% 2) Recovers independent k-space bands (asymmetric about centre) 
% 3) Completes k-space via conjugate symmetry
%
% The final image is then reconstructed by selecting the optimal off-
% resonance frequnecy in each voxel of the super-resolution image:
% 1) Low-pass filters k-space at each frequency offset
% 2) Uses 1) to generate smooth phase maps
% 3) Chooses the frequency offset which minimises the absolute phase
% 
% The final images reproduce the results shown in [REF].
%
% PJL 08/04/20

clear all

%% Load data & dependencies
addpath('./external/')             % Add path to ifft2c and fft2c functions
load('../data/complexFID.mat')     % Load demo dataset

N=size(complexFID);     % Determine the data dimensions
n_asym_bands=(N(3)/2)+1;    % Number of bands to recover by least squares
                            % *This code assumes image dimensions are even*

%% Decide on number of offsets for MF reconstruction
frq_offs=100;

%% Initialise arrays for S bands and MF recon
S_bands=zeros(N(1),N(2),n_asym_bands);
Im_mf=zeros(N(1)*N(3),N(2),frq_offs);

%% Loop through offsets for MF reconstruction
for frq=1:frq_offs

% Create array of pattern shift phases
phases=2*pi*((0:N(3)-1)/N(3))';
% Add an offset for MF reconstruction
phases=phases+frq/frq_offs*2*pi;

% Construct bands/shifts matrix psi 
psi=exp(-1i*(1:n_asym_bands).*phases);

% Obtain least squares solution, voxel by voxel
for pix_a=1:N(1)
    for pix_b=1:N(2)        
        S_bands(pix_a,pix_b,:)=psi\squeeze(complexFID(pix_a,pix_b,:)); 
    end
end

% Shift the bands to the correct location in k-space
for k=1:n_asym_bands
    S_asym((k-1)*N(1)+1:k*N(1),:)=S_bands(:,:,n_asym_bands+1-k);
end

% Complete S using conjugate symmetry
S = [S_asym(N(1)/2:end,:);rot90(conj(S_asym(N(1)/2+1:end-N(1)-1,:)),2)];

% Store the resulting image at this offset 
Im_mf(:,:,frq)=ifft2c(S);

end

%% Low-pass filter S at each offset to obtain smooth phase map
% Create mask for central band of S
centvox=N(1)*N(3)/2;
bandw=-N(1)/2:N(1)/2-1;
zf_mat=zeros(size(Im_mf(:,:,1)));
zf_mat(centvox+bandw,:)=1;

% Apply the mask to S at each offset and create low pass images
Im_zf=ifft2c(fft2c(Im_mf).*repmat(zf_mat,[1 1 frq_offs]));

% Determine optimal offset which minimises the absolute phase
[~,opt_off]=min(abs(angle(Im_zf)),[],3);

% Choose the optimal offset for each voxel in the image
Im_recon=zeros(size(Im_mf(:,:,1)));
for pix_a=1:size(Im_mf,1)
    for pix_b=1:size(Im_mf,2)        
        Im_recon(pix_a,pix_b)=Im_mf(pix_a,pix_b,mod(opt_off(pix_a,pix_b)-1,frq_offs)+1);
    end
end

%% Display the results
figure(1)
subplot(1,4,1)
imagesc(mean(abs(ifft2c(complexFID)),3));
caxis([0 1200])
set(gca,'DataAspectRatio',[1 1/6 1])
title({'Mean of','low-res. images'})
axis off
subplot(1,4,2)
imagesc(mean(abs(Im_zf),3));
caxis([0 250])
set(gca,'DataAspectRatio',[1 6 1])
title({'Mean of zerofilled','low-res. images'})
axis off
subplot(1,4,3)
imagesc(abs(Im_mf(:,:,1)));
caxis([0 250])
set(gca,'DataAspectRatio',[1 6 1])
title({'Super-resolution image','without MF recon'})
axis off
subplot(1,4,4)
imagesc(abs(Im_recon));
caxis([0 250])
set(gca,'DataAspectRatio',[1 6 1])
title({'Super-resolution image','with MF recon'})
axis off

colormap gray
set(gcf,'color','w')