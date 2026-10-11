function C = audit_config()
% Paths are relative to this public repository. Configure NIfTI paths only if rebuilding maps.
project = fileparts(fileparts(mfilename('fullpath')));
C.root = fullfile(project, 'results');
C.phenotypeInputDir = fullfile(project, 'data');
C.nmTargetFile = fullfile(project, 'data', 'nm_age6_80_75.mat');
C.dvTargetFile = fullfile(project, 'data', 'dv_1_age7_80_74.mat');
C.diseaseTargetFile = fullfile(project, 'data', 'disease_effect_sizes.mat');
C.networkFile = fullfile(project, 'data', 'molecular_nets_remove_dice_auto.mat');
C.spatialControlFile = fullfile(project, 'data', 'spa_dis.mat');
% Set these paths when rebuilding biological networks from NIfTI.
C.spm = '/path/to/spm12';
C.template = '/path/to/Neuromark_fMRI_1.0.nii';
C.ctDir = '/path/to/celltype';
C.ntDir = '/path/to/neurotransmitter';
C.mcDir = '/path/to/mitochondria';

C.nPerm = 1000;
C.seed = 20261008;
C.zThreshold = 3;
end
