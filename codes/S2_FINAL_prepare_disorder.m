function S2_FINAL_prepare_disorder()
% Prepare fifteen disease effect-size networks for spatial analysis.
% Five patient one-sample Cohen dz, then ten patient-vs-patient Cohen d.
% Spearman partial spatial control, node permutations,
C=audit_config();addpath(fileparts(mfilename('fullpath')));
N=load_network_inputs(C);
folder=fullfile(C.root,'multiplecorr');
assert(isfolder(folder),'Expected existing multiplecorr dir: %s',folder);
P=load(C.diseaseTargetFile,'D');D=P.D;
assert(isequal(size(D.raw_targets),[1378 15]),'Expected 15 final disease targets');
staging=fullfile(folder,'_FINAL_15_STAGING');
if isfolder(staging),error('Old staging exists; inspect/remove manually before retrying: %s',staging);end
mkdir(staging);
rng(C.seed,'twister');
mods={'ct','nt','mc'};
expected={... 
 {'L2_3IT','L4IT','L6IT','L6ITCar3','L5IT','L5ET','L5_6NP','L6CT','L6b','Chandelier','Pvalb','Sst','SstChodl','Lamp5','Sncg','Pax6','Lamp5Lhx6','Vip','Astro','Endo','OPC','VLMC','Micro_PVM','Oligo'}, ...
 {'5HT1a','5HT1b','5HT2a','5HT4','5HT6','5HTT','D1','D2','DAT','NET','H3','A4B2','M1','VAChT','CB1','MOR','NMDA','mGluR5','GABAa'}, ...
 {'CI','CII','CIV','MitoD','MRC','TRC'}};
for im=1:3
    md=mods{im};
    givenNames=N.names.(md);
    assert(isequal(cellstr(string(givenNames(:))),expected{im}(:)), ...
        'Unexpected marker order: %s',md);
    switch md
        case 'ct',nodes=6:49; Xnet=N.celltype_net(nodes,nodes,:);Z1=N.celltype_dice;
        case 'nt',nodes=1:53;Xnet=N.neurotransmitter_net;Z1=N.dice;
        case 'mc',nodes=1:53;Xnet=N.mitochondral_net;Z1=N.dice;
    end
    if strcmp(md,'ct'),Z2=N.celltype_spa_dis;else,Z2=N.spa_dis;end
        Z2=Z2-diag(diag(Z2));
    nCases=15;nFeatures=size(Xnet,3);
    rval=zeros(nCases,nFeatures);pval=zeros(nCases,nFeatures);
    quantile=zeros(nCases,nFeatures,6);
    for k=1:nCases
        Yfull=triangular_conversion(D.raw_targets(:,k),'vec2mat','lower');
        Ymat=Yfull(nodes,nodes);
        [x,y,yperm]=partial_mantel(Xnet,Ymat,Z1,Z2,C.nPerm);
        [rval(k,:),pval(k,:),quantile(k,:,:)]=single_corr(x,y,yperm);
        save(fullfile(staging,sprintf('disorder_%s_%d.mat',md,k)),'x','y','yperm');
        fprintf('  disorder / %s %2d/%d [%s]\n',md,k,nCases,D.target_names{k});
    end
    % Match S4 multivariable families: Fig4 (5 x markers) and Fig5 (10 x markers).
    fdr=nan(size(pval));
    p4=pval(1:5,:);p5=pval(6:15,:);
    fdr(1:5,:)=reshape(mafdr(p4(:),'BHFDR',true),size(p4));
    fdr(6:15,:)=reshape(mafdr(p5(:),'BHFDR',true),size(p5));
    names=N.names.(md);
    meta=struct('task','disorder','modality',md,'names',{names},'nCases',15, ...
         'nPerm',C.nPerm,'targetNames',{D.target_names},'targetKinds',{D.target_kinds},...
         'fdr_family','separate_by_modality; Fig4 rows1to5; Fig5 rows6to15');
    save(fullfile(staging,sprintf('SINGLEVAR_disorder_%s.mat',md)), ...
         'rval','pval','fdr','quantile','meta');
    T=table((1:nFeatures)',string(names(:)),'VariableNames',{'column','marker'});
    writetable(T,fullfile(staging,sprintf('COLUMN_ORDER_disorder_%s.csv',md)));
end
save(fullfile(staging,'FINAL_DISEASE_NETWORK_DATA.mat'),'D');
manifest=table({'disorder';'disorder';'disorder'},mods(:),repmat(15,3,1), ...
    'VariableNames',{'task','modality','cases'});
writetable(manifest,fullfile(staging,'RUN_MANIFEST.csv'));
backup=fullfile(folder,['_SUPERSEDED_DISORDER_' datestr(now,'yyyymmdd_HHMMSS')]);
mkdir(backup);
patterns={'disorder_*.mat','results_disorder_*.mat', ...
 'SINGLEVAR_disorder_*.mat','COLUMN_ORDER_disorder_*.csv', ...
 'plot_results_disorder_*.mat','plot_results_dcomparision_*.mat', ...
 'NEW_MARKERS_disorder_*.csv','NEW_MARKERS_dcomparision_*.csv', ...
 'FINAL_DISEASE_NETWORK_DATA.mat','RUN_MANIFEST.csv'};
for z=1:numel(patterns)
    files=dir(fullfile(folder,patterns{z}));
    for q=1:numel(files)
        if files(q).isdir,continue;end
        [ok,msg]=movefile(fullfile(files(q).folder,files(q).name),fullfile(backup,files(q).name));
        assert(ok,'Unable to archive previous disease file: %s',msg);
    end
end
files=dir(fullfile(staging,'*'));
for q=1:numel(files)
    if files(q).isdir,continue;end
    [ok,msg]=movefile(fullfile(files(q).folder,files(q).name),fullfile(folder,files(q).name));
    assert(ok,'Unable to activate final disease input: %s',msg);
end
rmdir(staging);
fid=fopen(fullfile(folder,'FINAL_15_STATUS.txt'),'w');
fprintf(fid,'S2_DONE: 15 disease targets; R and S4 NOT YET DONE\n');
fclose(fid);
fprintf('\nS2 COMPLETE: 15 NEW disease targets x 3 modalities.\n');
fprintf('Old disease files archived at %s\n',backup);
fprintf('Existing nm and dv_1 files were NOT modified.\n');
end

function [x_resid,y_resid,y_perm]=partial_mantel(X,Y,Z1,Z2,nPerm)
% squareform() uses identical edge ordering for x,y,Z and permuted Y.
num=size(X,3); e=size(X,1)*(size(X,1)-1)/2;
x=zeros(num,e);
for q=1:num
    mat=X(:,:,q);x(q,:)=squareform(mat);
end
y=squareform(Y);
z1=squareform(Z1);z2=squareform(Z2);
Z=[z1;z2];
vecA_resid=partial_remove_covs(x',y',Z','type','spearman');
YY=triangular_conversion(vecA_resid(:,end),'vec2mat','lower');
n=size(X,1);
y_perm=zeros(nPerm,e);
for q=1:nPerm
    p=randperm(n);
    Yp=YY(p,p);
    y_perm(q,:)=squareform(Yp);
end
x_resid=vecA_resid(:,1:end-1);
y_resid=vecA_resid(:,end);
y_perm=y_perm';
end
function [r,p,qt]=single_corr(X,Y,perms)
nPerm=size(perms,2);
r=corr(X,Y,'type','spearman');
rPerm=corr(X,perms,'type','spearman');
p=(sum(abs(rPerm)>=abs(r),2))'/nPerm;
r=r';
percentiles=[99 95 75 25 5 1];
qt=zeros(1,size(X,2),numel(percentiles));
for q=1:size(X,2)
    qt(1,q,:)=prctile(rPerm(q,:),percentiles);
end
end
