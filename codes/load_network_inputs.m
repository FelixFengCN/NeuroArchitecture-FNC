function N = load_network_inputs(C)
% Load released similarity networks and spatial covariates directly.
N = load(C.networkFile, 'celltype_net','neurotransmitter_net','mitochondral_net');
S = load(C.spatialControlFile,'celltype_dice','celltype_spa_dis','dice','spa_dis');
assert(isequal(size(N.celltype_net),[53 53 24]));
assert(isequal(size(N.neurotransmitter_net),[53 53 19]));
assert(isequal(size(N.mitochondral_net),[53 53 6]));
assert(isequal(size(S.celltype_dice),[44 44]));
assert(isequal(size(S.celltype_spa_dis),[44 44]));
assert(isequal(size(S.dice),[53 53]));
assert(isequal(size(S.spa_dis),[53 53]));
N.celltype_dice=S.celltype_dice;
N.celltype_spa_dis=S.celltype_spa_dis;
N.dice=S.dice;
N.spa_dis=S.spa_dis;
N.names.ct={'L2_3IT','L4IT','L6IT','L6ITCar3','L5IT','L5ET','L5_6NP','L6CT','L6b','Chandelier','Pvalb','Sst','SstChodl','Lamp5','Sncg','Pax6','Lamp5Lhx6','Vip','Astro','Endo','OPC','VLMC','Micro_PVM','Oligo'};
N.names.nt={'5HT1a','5HT1b','5HT2a','5HT4','5HT6','5HTT','D1','D2','DAT','NET','H3','A4B2','M1','VAChT','CB1','MOR','NMDA','mGluR5','GABAa'};
N.names.mc={'CI','CII','CIV','MitoD','MRC','TRC'};
end
