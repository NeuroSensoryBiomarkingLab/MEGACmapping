%script to apply a t-test to the beta files from spm_vol
%must be placed in same directory as beta files
%authored by Kuwook Cha at the Montreal Neurological Institute for use by the EVS lab

clear Conditions

%create paths to all beta files 
%store in niiFilePath field of Conditions structure
for ic = 1:8
    Conditions(ic).niiFilePath{1} = sprintf('beta_%04d.nii',ic);
    Conditions(ic).niiFilePath{2} = sprintf('beta_%04d.nii',ic+14);
    Conditions(ic).niiFilePath{3} = sprintf('beta_%04d.nii',ic+28);
end

%extract beta file matrices into data field of Conditions structure
for ic = 1:length(Conditions)
    Conditions(ic).data(:,:,:,1) = niftiread(Conditions(ic).niiFilePath{1});
    Conditions(ic).data(:,:,:,2) = niftiread(Conditions(ic).niiFilePath{2});
    Conditions(ic).data(:,:,:,3) = niftiread(Conditions(ic).niiFilePath{3});
end

%% pseudo-t-test

allBetas = cat(5,Conditions(:).data);

%store matrix from ResMs.nii into ResMs
ResMs = niftiread('ResMs.nii');

%take the mean of all beta values across the 4th dimmension 
nom = squeeze(mean(allBetas,4));

%create a matrix of t statistics
t = nom./sqrt(ResMs);

df = 2*8*2-8; %degrees of freedom
p = 1-tcdf(t,df);
h = any((p<1 & t>1.4),4); % adjust p value here

%display data
figure;

isl = 38; % slice (z)
% subplot(2,2,1)
% imagesc(t(:,:,isl));
% subplot(2,2,2)
% imagesc(-log10(p(:,:,isl)))
% subplot(2,2,3)
% imagesc(h(:,:,isl));

tun = [];
%extract beta values across frequencies
for ic = 1:8    
    beta = mean(allBetas(:,:,:,:,ic),4);
    tun(ic,:) = beta(h);
end

[~,bf] = max(tun); % fill bf with max of tuning function (best frequency)
figure;hist(bf,[1:8]);title('Subject Tuning Frequency');xlabel('Stimulus Number');ylabel('Number of Voxels');

%fill bfMap with preferred frequency (full NIfTI matrix)
bfMap = zeros(size(h));
bfMap(h==1) = bf;
figure;imagesc(bfMap(:,:,isl));colormap(jet)

% write bfMap to NIfTI file
betaHeader = spm_vol('beta_0001.nii'); %steal header from compatible NIfTI
betaHeader.fname = 'bfm_t14.nii'; %name for tuning map file
betaHeader.private.dat.fname = betaHeader.fname;
spm_write_vol(betaHeader,bfMap);





