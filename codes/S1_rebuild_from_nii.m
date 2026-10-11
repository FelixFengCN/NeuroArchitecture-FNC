function S1_rebuild_from_nii()
% Construct cell-type, neurotransmitter, and mitochondrial
% microprofile similarity networks from NIfTI atlases.

C = audit_config();

addpath(genpath(C.spm));
addpath(fileparts(mfilename('fullpath')));

assert(exist('spm_vol','file') ~= 0, ...
    'SPM12 is required.');

if ~isfolder(C.root)
    mkdir(C.root);
end

%% 1. NeuroMark ICN masks

V = spm_vol(C.template);

assert(numel(V) == 53, ...
    'The NeuroMark template must contain 53 ICNs.');

T = spm_read_vols(V);

assert(size(T,4) == 53, ...
    'Expected a 4D NeuroMark template.');

T(abs(T) < 3) = 0;

mask = T ~= 0;

refdim = size(T(:,:,:,1));
refaffine = V(1).mat;

%% 2. Biological marker names and order

ctRaw = { ...
    'L2_3IT','L4IT','L6IT','L6ITCar3','L5IT', ...
    'Lamp5','Sncg','Pax6','Lamp5Lhx6','Vip', ...
    'L5ET','L5_6NP','L6CT','L6b','Chandelier', ...
    'Pvalb','Sst','SstChodl','Astro','Endo', ...
    'OPC','VLMC','Micro_PVM','Oligo'};

ctIndex = [1:5,11:14,15:18,6:10,19:24];
ctNames = ctRaw(ctIndex);

ntNames = { ...
    '5HT1a','5HT1b','5HT2a','5HT4','5HT6','5HTT', ...
    'D1','D2','DAT','NET','H3','A4B2','M1', ...
    'VAChT','CB1','MOR','NMDA','mGluR5','GABAa'};

mcNames = {'CI','CII','CIV','MitoD','MRC','TRC'};

names = struct( ...
    'ct',{ctNames}, ...
    'nt',{ntNames}, ...
    'mc',{mcNames});

%% 3. Voxel sampling

% AHBA valid-voxel mask for cell-type maps
ahba = read_nii( ...
    fullfile(C.ctDir,'ahba_mask.nii'), ...
    refdim,refaffine);

ahbaMask = ahba ~= 0;

% Assign a unique label to each AHBA voxel
atlas_celltype = zeros(refdim);

atlas_celltype(ahbaMask) = 1:nnz(ahbaMask);

ctSample = cell(44,1);
ntSample = cell(53,1);

% CT: cortical ICNs 6–49
for ic = 6:49

    temp = atlas_celltype;
    temp(~mask(:,:,:,ic)) = 0;

    [~, indices, ~] = labelled_voxel_centroids(temp);

    ctSample{ic-5} = indices;
end

% NT and MC: all 53 ICNs
for ic = 1:53

    temp = T(:,:,:,ic);
    nonzero = temp ~= 0;

    temp(nonzero) = 1:nnz(nonzero);

    [~, indices, ~] = labelled_voxel_centroids(temp);

    ntSample{ic} = indices;
end

%% 4. Cell-type networks

celltype_net = zeros(53,53,24);

for k = 1:24

    filepath = fullfile(C.ctDir,[ctNames{k} '.nii']);

    A = read_nii(filepath,refdim,refaffine);
    A = A(:);

    samples = cellfun( ...
        @(ind) A(ind),ctSample, ...
        'UniformOutput',false);

    S = morphometric_similarity( ...
        samples,[],0,'one_minus_sqrt_jsd');

    celltype_net(6:49,6:49,k) = S-eye(44);

    fprintf('CT %2d/24: %s\n',k,ctNames{k});
end

%% 5. Neurotransmitter networks

neurotransmitter_net = zeros(53,53,19);

for k = 1:19

    filepath = fullfile(C.ntDir,[ntNames{k} '.nii']);

    A = read_nii(filepath,refdim,refaffine);
    A(~isfinite(A)) = 0;
    A = A(:);

    samples = cellfun( ...
        @(ind) A(ind),ntSample, ...
        'UniformOutput',false);

    S = morphometric_similarity( ...
        samples,[],0,'one_minus_sqrt_jsd');

    neurotransmitter_net(:,:,k) = S-eye(53);

    fprintf('NT %2d/19: %s\n',k,ntNames{k});
end

%% 6. Mitochondrial networks

mitochondral_net = zeros(53,53,6);

for k = 1:6

    filepath = fullfile(C.mcDir,[mcNames{k} '.nii.gz']);

    A = read_nii(filepath,refdim,refaffine);
    A(isnan(A)) = 0;
    A = A(:);

    samples = cellfun( ...
        @(ind) A(ind),ntSample, ...
        'UniformOutput',false);

    S = morphometric_similarity( ...
        samples,[],0,'one_minus_sqrt_jsd');

    mitochondral_net(:,:,k) = S-eye(53);

    fprintf('MC %d/6: %s\n',k,mcNames{k});
end

%% 7. Spatial covariates

% Compute Dice overlap using the actual sampled voxels
celltype_dice = make_dice(ctSample);
dice = make_dice(ntSample);

% Load spatial distance matrices
L = load(C.spatialControlFile, ...
    'spa_dis','celltype_spa_dis');

spa_dis = double(L.spa_dis);
celltype_spa_dis = double(L.celltype_spa_dis);

assert(isequal(size(celltype_dice),[44 44]));
assert(isequal(size(dice),[53 53]));
assert(isequal(size(celltype_spa_dis),[44 44]));
assert(isequal(size(spa_dis),[53 53]));

celltype_spa_dis(1:45:end) = 0;
spa_dis(1:54:end) = 0;

%% 8. Save networks

save(fullfile(C.root,'nets_from_nii.mat'), ...
    'celltype_net', ...
    'neurotransmitter_net', ...
    'mitochondral_net', ...
    'celltype_dice', ...
    'celltype_spa_dis', ...
    'dice', ...
    'spa_dis', ...
    'names', ...
    'ctIndex');

fprintf('\nNetwork construction completed.\n');
fprintf('Saved: %s\n', ...
    fullfile(C.root,'nets_from_nii.mat'));

end


%% Helper: extract labelled voxel centroids

function [xyz,idx,labels] = labelled_voxel_centroids(label_volume)

labels = unique(label_volume(:));
labels(labels == 0) = [];

idx = find(label_volume ~= 0);

[assigned_labels,order] = sort(label_volume(idx));
idx = idx(order);

assert(isequal(double(assigned_labels(:)),double(labels(:))), ...
    'Each nonzero voxel must have a unique label.');

[x,y,z] = ind2sub(size(label_volume),idx);

xyz = [x(:),y(:),z(:)];

idx = sub2ind(size(label_volume), ...
    xyz(:,1),xyz(:,2),xyz(:,3));

end


%% Helper: read NIfTI

function A = read_nii(filepath,refdim,refaffine)

assert(isfile(filepath), ...
    'Missing NIfTI file: %s',filepath);

if endsWith(lower(filepath),'.nii.gz')

    tmp = tempname;
    mkdir(tmp);

    cleanupObj = onCleanup(@() rmdir(tmp,'s'));

    files = gunzip(filepath,tmp);
    filepath = files{1};
end

V = spm_vol(filepath);

assert(numel(V) == 1, ...
    'Expected a 3D NIfTI file.');

assert(isequal(double(V.dim(1:3)),double(refdim)), ...
    'NIfTI dimensions do not match NeuroMark.');

assert(max(abs(V.mat(:)-refaffine(:))) < 1e-3, ...
    'NIfTI affine does not match NeuroMark.');

A = spm_read_vols(V);

end


%% Helper: Dice overlap

function dice = make_dice(samples)

N = numel(samples);
dice = zeros(N,N);

for i = 1:N

    a = unique(samples{i});

    for j = i+1:N

        b = unique(samples{j});

        overlap = 2*numel(intersect(a,b));
        total = numel(a)+numel(b);

        if total > 0

            dice(i,j) = overlap/total;
            dice(j,i) = dice(i,j);

        end
    end
end

end
