function S2_FINAL_prepare_age_from_mat()

C=audit_config();addpath(fileparts(mfilename('fullpath')));
N=load_network_inputs(C);
folder=fullfile(C.root,'multiplecorr');if ~isfolder(folder),mkdir(folder);end
mods={'ct','nt','mc'};
A=load(C.nmTargetFile,'fnc_nm','ages');
B=load(C.dvTargetFile,'fnc_dv_1','ages');
assert(isequal(size(A.fnc_nm),[1378 75]) && isequal(A.ages(:),(6:80)'), 'NM MAT mismatch');
assert(isequal(size(B.fnc_dv_1),[1378 74]) && isequal(B.ages(:),(7:80)'), 'DV MAT mismatch');
for it=1:2
    if it==1,task='nm';Yall=A.fnc_nm;else,task='dv_1';Yall=B.fnc_dv_1;end
    nCases=size(Yall,2);
    for im=1:3
        md=mods{im};
        switch md
            case 'ct',nodes=6:49; Xnet=N.celltype_net(nodes,nodes,:);Z1=N.celltype_dice;
            case 'nt',nodes=1:53; Xnet=N.neurotransmitter_net;Z1=N.dice;
            case 'mc',nodes=1:53; Xnet=N.mitochondral_net;Z1=N.dice;
        end
        if strcmp(md,'ct'),Z2=N.celltype_spa_dis;else,Z2=N.spa_dis;end
        Z2=Z2-diag(diag(Z2));
        names=N.names.(md);nFeatures=size(Xnet,3);
        rval=zeros(nCases,nFeatures);pval=rval;quantile=zeros(nCases,nFeatures,6);
        rng(C.seed,'twister');
        for k=1:nCases
            Yfull=triangular_conversion(Yall(:,k),'vec2mat','lower');
            Ymat=Yfull(nodes,nodes);
            [x,y,yperm]=partial_mantel(Xnet,Ymat,Z1,Z2,C.nPerm);
            [rval(k,:),pval(k,:),quantile(k,:,:)]=single_corr(x,y,yperm);
            save(fullfile(folder,sprintf('%s_%s_%d.mat',task,md,k)),'x','y','yperm');
            fprintf('%s/%s %d/%d\n',task,md,k,nCases);
        end
        fdr=reshape(mafdr(pval(:),'BHFDR',true),size(pval));
        meta=struct('task',task,'modality',md,'names',{names},'nCases',nCases,'nPerm',C.nPerm);
        save(fullfile(folder,sprintf('SINGLEVAR_%s_%s.mat',task,md)), ...
             'rval','pval','fdr','quantile','meta');
    end
end
fprintf('Age S2 complete (NM 75; DV1 74). Run R age stage next.\n');
end

function [x_resid,y_resid,y_perm]=partial_mantel(X,Y,Z1,Z2,nPerm)
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
