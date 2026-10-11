function S4_FINAL_postprocess()

C=audit_config();addpath(fileparts(mfilename('fullpath')));folder=fullfile(C.root,'multiplecorr');
N=load_network_inputs(C);
P=load(fullfile(folder,'FINAL_DISEASE_NETWORK_DATA.mat'),'D');D=P.D;
assert(isequal(size(D.raw_targets),[1378 15]),'Wrong disease targets');
mods={'ct','nt','mc'};
for im=1:3
    md=mods{im};names=N.names.(md);nFeat=numel(names);
    for task={'nm','dv_1'}
        taskName=task{1};
        if strcmp(taskName,'nm'),ks=1:75;else,ks=1:74;end
        S=read_R_batch(folder,taskName,md,ks,C.nPerm);
        assert(size(S.weis,2)==nFeat,'Wrong marker count %s %s',taskName,md);
        S.marker_names=names;S.target_names=arrayfun(@(i)sprintf('%s_%02d',taskName,i),ks,'UniformOutput',false);
        save(fullfile(folder,sprintf('plot_results_%s_%s.mat',taskName,md)),'-struct','S');
        emit_csv(S,names,taskName,md,S.target_names,folder);
        fprintf('Restored saved-R summary %s/%s (%d age points)\n',taskName,md,numel(ks));
    end

    for k=1:15
        I=load(fullfile(folder,sprintf('disorder_%s_%d.mat',md,k)),'y');
        O=load(fullfile(folder,sprintf('results_disorder_%s_%d.mat',md,k)),'observed_y');
        assert(numel(I.y)==numel(O.observed_y) && isequaln(I.y(:),O.observed_y(:)), ...
            'R observed_y differs from final MATLAB target in %s, #%d',md,k);
    end
    Sfirst=read_R_batch(folder,'disorder',md,1:5,C.nPerm);
    Ssecond=read_R_batch(folder,'disorder',md,6:15,C.nPerm);
    assert(size(Sfirst.weis,2)==nFeat && size(Ssecond.weis,2)==nFeat, 'Wrong predictor count');
    Sfirst.marker_names=names;Sfirst.target_names=D.target_names(1:5);
    Sfirst.phenotype='One-sample Cohen dz = mean patient dscore / sample SD';
    Ssecond.marker_names=names;Ssecond.target_names=D.target_names(6:15);
    Ssecond.phenotype='Between-patient-group pooled-SD Cohen d, first minus second';
    save(fullfile(folder,sprintf('plot_results_disorder_%s.mat',md)),'-struct','Sfirst');
    save(fullfile(folder,sprintf('plot_results_dcomparision_%s.mat',md)),'-struct','Ssecond');
    emit_csv(Sfirst,names,'disorder',md,D.target_names(1:5),folder);
    emit_csv(Ssecond,names,'dcomparision',md,D.target_names(6:15),folder);
    fprintf('Completed FINAL 5 dz + 10 pairwise d for %s\n',md);
end
S5_FINAL_update_disorder_singlevar_FDR(); % correct disease univariate FDR too
fid=fopen(fullfile(folder,'FINAL_15_STATUS.txt'),'w');
fprintf(fid,'S4_DONE: 15 disease targets; R verified; disease singlevar and multivar use Fig4/Fig5 separate BH families\n');fclose(fid);
fprintf('FINAL S4 COMPLETE: 12 flat plotting MAT files; no historical outputs used.\n');
end

function S=read_R_batch(folder,task,md,ks,nPerm)
% Exact previous A09/A10 aggregation, with original > (not >= or +1).
n=numel(ks);S=struct;
for t=1:n
    k=ks(t);f=fullfile(folder,sprintf('results_%s_%s_%d.mat',task,md,k));
    assert(isfile(f),'Missing R output: %s',f);
    R=load(f,'adj_r2','weis');adj=R.adj_r2(:);w=R.weis;
    if size(w,1)~=nPerm+1 && size(w,2)==nPerm+1,w=w';end
    assert(numel(adj)==nPerm+1 && size(w,1)==nPerm+1,'Bad permutation dimensions in %s',f);
    nf=size(w,2);
    if t==1
        S.r2=zeros(n,1);S.r2_null=zeros(n,nPerm);S.r2_p=zeros(n,1);
        S.weis=zeros(n,nf);S.weis_null=zeros(n,nPerm,nf);S.weis_p=zeros(n,nf);
        S.sub_r2=zeros(n,nf);S.sub_r2_null=zeros(n,nPerm,nf);S.sub_r2_p=zeros(n,nf);
    end
    assert(size(S.weis,2)==nf,'Predictor count changed: %s',f);
    S.r2(t)=adj(1);S.r2_null(t,:)=adj(2:end)';
    S.r2_p(t)=sum(S.r2_null(t,:)>S.r2(t))/nPerm;
    S.weis(t,:)=w(1,:);S.weis_null(t,:,:)=w(2:end,:);
    S.weis_p(t,:)=sum(w(2:end,:)>w(1,:),1)/nPerm;
    S.sub_r2(t,:)=w(1,:)*adj(1);
    nullContrib=w(2:end,:).*adj(2:end);
    S.sub_r2_null(t,:,:)=nullContrib;
    S.sub_r2_p(t,:)=sum(nullContrib>S.sub_r2(t,:),1)/nPerm;
end
S.r2_fdr=mafdr(S.r2_p,'BHFDR',true);
S.sub_r2_fdr=reshape(mafdr(S.sub_r2_p(:),'BHFDR',true),size(S.sub_r2_p));
S.weis_zscore=zscore(S.weis');
end

function emit_csv(S,names,task,md,targetNames,folder)
n=numel(targetNames);p=numel(names);
T=table('Size',[n*p,10], ...
 'VariableTypes',{'string','string','string','double','string','double','double','double','double','double'}, ...
 'VariableNames',{'task','modality','target','marker_column','marker','adjusted_R2','R2_p','R2_FDR','genizi_weight','marker_FDR'});
q=0;
for k=1:n
    for j=1:p
        q=q+1;T.task(q)=string(task);T.modality(q)=string(md);
        T.target(q)=string(targetNames{k});T.marker_column(q)=j;
        T.marker(q)=string(names{j});T.adjusted_R2(q)=S.r2(k);
        T.R2_p(q)=S.r2_p(k);T.R2_FDR(q)=S.r2_fdr(k);
        T.genizi_weight(q)=S.weis(k,j);T.marker_FDR(q)=S.sub_r2_fdr(k,j);
    end
end
writetable(T,fullfile(folder,sprintf('FINAL_MARKERS_%s_%s.csv',task,md)));
end
