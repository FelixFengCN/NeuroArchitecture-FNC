function S5_FINAL_update_disorder_singlevar_FDR()
% Recompute only disease SINGLE-VARIATE BH q-values from SAVED permutation p-values.
% Match the multivariable families used in S4: Fig4 (5 x P), Fig5 (10 x P).
% CT/NT/MC are separate correction families; NO permutations, S2, R, or S4 rerun.
C=audit_config();folder=fullfile(C.root,'multiplecorr');
assert(isfolder(folder),'Missing multiplecorr directory: %s',folder);
A=load(fullfile(folder,'FINAL_DISEASE_NETWORK_DATA.mat'),'D');D=A.D;
assert(isequal(size(D.raw_targets),[1378,15]),'Missing final 15 Cohen-d targets');
mods={'ct','nt','mc'};nMarkers=[24,19,6];
for im=1:3
    md=mods{im};p=nMarkers(im);
    f=fullfile(folder,sprintf('SINGLEVAR_disorder_%s.mat',md));
    assert(isfile(f),'Missing disease single-variable results: %s',f);
    S=load(f);
    assert(isfield(S,'pval') && isfield(S,'fdr') && isequal(size(S.pval),[15 p]), ...
        'Invalid single-variable p matrix: %s',f);
    assert(all(isfinite(S.pval(:))) && all(S.pval(:)>=0 & S.pval(:)<=1), ...
        'Non-finite/out-of-range permutation p-values: %s',f);
    assert(isfield(S,'meta') && isfield(S.meta,'targetNames') && ...
        isequal(cellstr(string(S.meta.targetNames(:))),cellstr(string(D.target_names(:)))), ...
        'Disease target names do not match FINAL 5+10 targets in %s',f);
    old_fdr=S.fdr;
    p4=S.pval(1:5,:);p5=S.pval(6:15,:);
    fdr=nan(size(S.pval));
    fdr(1:5,:)=reshape(mafdr(p4(:),'BHFDR',true),size(p4));
    fdr(6:15,:)=reshape(mafdr(p5(:),'BHFDR',true),size(p5));
    S.fdr=fdr;
    S.meta.fdr_family='separate_by_modality; Fig4 rows1to5; Fig5 rows6to15';
    assert(isfield(S,'rval') && isequal(size(S.rval),size(fdr)), ...
        'Missing or unexpected rval shape: %s',f);
    if isfield(S,'quantile'),assert(size(S.quantile,1)==15 && size(S.quantile,2)==p);end
    assert(all(isfinite(fdr(:))) && all(fdr(:)>=0 & fdr(:)<=1));
    save(f,'-struct','S'); % MATLAB default MAT version, preserve all other fields
    fprintf('FDR FIX %s: Fig4 %d tests, Fig5 %d tests; changed q=%d; stars Fig4 %d->%d, Fig5 %d->%d\n', ...
        md,5*p,10*p,nnz(abs(old_fdr(:)-fdr(:))>1e-12), ...
        nnz(old_fdr(1:5,:)<0.05),nnz(fdr(1:5,:)<0.05), ...
        nnz(old_fdr(6:15,:)<0.05),nnz(fdr(6:15,:)<0.05));
end
fprintf('All disease univariate q-values now use EXACT S4 figure/modality BH families.\n');
end
