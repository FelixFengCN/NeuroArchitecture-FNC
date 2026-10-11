function D = FINAL_build_disease_phenotypes(C)
% Five one-sample Cohen's dz and ten between-patient-group Cohen's d.
% All analyses use patient D-scores only.

if nargin == 0
    C = audit_config();
end

A = load(C.individualDscoreFile, 'dscore', 'covs_disorder');

assert(isfield(A, 'dscore') && isfield(A, 'covs_disorder'), ...
    'Missing individual D-score or covariates');

Y = double(A.dscore);
Z = double(A.covs_disorder);

assert(size(Y,2) == 1378 && ...
       size(Y,1) == size(Z,1) && ...
       size(Z,2) >= 4, ...
       'Individual D-score dimension mismatch');

assert(all(isfinite(Y(:))), ...
    'NaN/Inf D-score values');

site = Z(:,3);
diagnosis = Z(:,4);

names = {'SCZ','BP_LA5C','MDD_CMDD','ASD','ADHD'};

siteSets = {6:13, 14, 15:18, 23:58, 59:60};

expectedN = [219 46 256 764 230];

patients = cell(1,5);
pn = zeros(1,5);

dz = zeros(1378,5);
t1 = zeros(1378,5);
p1 = zeros(1378,5);
q1 = zeros(1378,5);

meanD = zeros(1378,5);
sdD = zeros(1378,5);

%% 1. One-sample patient Cohen's dz

for g = 1:5

    in = ismember(site,siteSets{g}) & ...
         ~ismember(site,[32 38 65]) & ...
         diagnosis == 1;

    patients{g} = Y(in,:);

    pn(g) = size(patients{g},1);

    assert(pn(g) == expectedN(g), ...
        '%s patient count mismatch: %d (expected %d)', ...
        names{g},pn(g),expectedN(g));

    meanD(:,g) = mean(patients{g},1)';
    sdD(:,g) = std(patients{g},0,1)';

    assert(all(isfinite(sdD(:,g)) & sdD(:,g)>0), ...
        'Zero/invalid patient SD for %s',names{g});

    dz(:,g) = meanD(:,g)./sdD(:,g);

    [~,p,~,stats] = ttest(patients{g},0);

    p1(:,g) = p(:);
    t1(:,g) = stats.tstat(:);

    q1(:,g) = mafdr(p(:),'BHFDR',true);

    assert(max(abs(t1(:,g)-dz(:,g)*sqrt(pn(g))))<1e-8, ...
        'Cohen dz / one-sample t mismatch');

    fprintf('ONE-SAMPLE %-9s n=%d dz on 1378 edges, q<.05=%d\n', ...
        names{g},pn(g),nnz(q1(:,g)<0.05));
end

%% 2. Between-patient-group Cohen's d

pairIndex = nchoosek(1:5,2);

pairNames = cell(1,10);

pairD = zeros(1378,10);
pairT = zeros(1378,10);
pairP = zeros(1378,10);
pairQ = zeros(1378,10);

for k = 1:10

    i = pairIndex(k,1);
    j = pairIndex(k,2);

    X = patients{i};
    W = patients{j};

    n1 = size(X,1);
    n2 = size(W,1);

    pooled = sqrt( ...
        ((n1-1)*var(X,0,1) + ...
         (n2-1)*var(W,0,1)) / ...
        (n1+n2-2));

    assert(all(isfinite(pooled) & pooled>0), ...
        'Invalid pooled SD in pair %d',k);

    pairD(:,k) = ...
        ((mean(X,1)-mean(W,1))./pooled)';

    [~,p,~,s] = ttest2(X,W);

    pairT(:,k) = s.tstat(:);
    pairP(:,k) = p(:);

    pairQ(:,k) = mafdr(p(:),'BHFDR',true);

    pairNames{k} = sprintf('%s_vs_%s', ...
        names{i},names{j});

    fprintf('PAIR %2d %-23s q<.05=%d\n', ...
        k,pairNames{k},nnz(pairQ(:,k)<0.05));
end

%% 3. Assemble 15 disease targets

rawTargets = [dz pairD];

targetNames = ...
    [strcat(names,'_one_sample_dz') pairNames];

targetKinds = ...
    [repmat({'ONE_SAMPLE_COHEN_DZ'},1,5) ...
     repmat({'BETWEEN_DISEASE_COHEN_D'},1,10)];

D = struct();

D.group_names = names;
D.n_patient = pn;

D.one_sample_dz = dz;
D.one_sample_t = t1;
D.one_sample_p = p1;
D.one_sample_q = q1;

D.patient_mean_dscore = meanD;
D.patient_sd_dscore = sdD;

D.pair_indices = pairIndex;
D.pair_names = pairNames;

D.pair_cohen_d = pairD;
D.pair_t = pairT;
D.pair_p = pairP;
D.pair_q = pairQ;

D.raw_targets = rawTargets;
D.target_names = targetNames;
D.target_kinds = targetKinds;

%% 4. Network matrices

D.one_sample_dz_net = to_network(dz);
D.one_sample_t_net = to_network(t1);
D.one_sample_p_net = to_network(p1);
D.one_sample_q_net = to_network(q1);
D.one_sample_sig_net = to_network(double(q1<0.05));

D.pair_cohen_d_net = to_network(pairD);
D.pair_t_net = to_network(pairT);
D.pair_p_net = to_network(pairP);
D.pair_q_net = to_network(pairQ);
D.pair_sig_net = to_network(double(pairQ<0.05));

D.bh_family = ...
    '1378 edges separately within each of the 5 or 10 contrasts';

D.effect_direction = ...
    'one-sample patient mean/SD; pair first-minus-second pooled-SD';

end

function nets = to_network(X)

nets = zeros(53,53,size(X,2));

for i = 1:size(X,2)
    nets(:,:,i) = triangular_conversion( ...
        X(:,i),'vec2mat','lower');
end

end