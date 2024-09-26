clear

% Load the dataset
parameters=load('/media/koba/MULTIBOOT/blindness_gradients/source/parameters.mat');
parameters=parameters.parameters;



% gradients_all=zeros(size(parameters,2), size(parameters(1).gradients,1),size(parameters(1).gradients,2));
for i=1:size(parameters,2)
    % gradients_sp80_na(i,:,:)=parameters(i).gradients_sp80_na;
    % gradients_sp90_na(i,:,:)=parameters(i).gradients_sp90_na;
    gradients_sp90_cs(i,:,:)=parameters(i).gradients_sp90_cs;
end
gradients_all=gradients_sp90_cs;
% Remove the outliers and check inter-intra-class correlations
cors=corr(gradients_all(:,:,1)');
low_cors=(abs(cors)<0.4);
low_cors_id=(sum(low_cors)>20)';
nansub=(sum(isnan(cors))>1)';

% gradients_all(logical(low_cors_id),:,:)=[];
parameters(logical(low_cors_id +nansub))=[];

% Load parcellation-related info
addpath(genpath('/home/koba/freesurfer/matlab'))
[~, label_l, colortable_l]=read_annotation('/home/koba/micapipe/parcellations/lh.glasser-360_mics.annot');
[~, label_r, colortable_r]=read_annotation('/home/koba/micapipe/parcellations/rh.glasser-360_mics.annot');
for i=1:size(colortable_l.table,1)
    label_l(label_l==colortable_l.table(i,5))=i;
    label_r(label_r==colortable_r.table(i,5))=i;
end
label_r=label_r+180;
label_r(label_r==181)=1;
glasser_label=[label_l;label_r]-1;
rmpath(genpath('/home/koba/freesurfer/matlab'))
% glasser_label=table2array(readtable('/media/koba/MULTIBOOT/blindness_gradients/source/parcellations/glasser_fsaverage5_labels.txt'));
glasser_metadata=readtable('/media/koba/MULTIBOOT/blindness_gradients/source/parcellations/HCP-MMP1_UniqueRegionList.csv');
glassercontinouscolors=table2array(readtable('/media/koba/MULTIBOOT/blindness_gradients/source/parcellations/glasser_continous_colors.txt'));
glassercontinouscolors=glassercontinouscolors(:,1:3);
addpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/BrainStat-master'))
reference_gradients=fetch_gradients();
[surf_lh, surf_rh] = fetch_template_surface('fsaverage5','layer','inflated');
[surf_lh_pial, surf_rh_pial] = fetch_template_surface('fsaverage5','layer','pial');
[network_names, yeo_colormap] = fetch_yeo_networks_metadata(7);
rmpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/BrainStat-master'))
addpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/ENIGMA-2.0.0'))


network_names={'visual', 'somatosensory','dorsal attention','ventral attention','limbic','frontoparietal','default mode'};
yeo_in_glasser=zeros(360,1);
for i=1:7
    yeo_in_glasser=yeo_in_glasser+strcmpi(glasser_metadata.community_yeo,network_names(i))*i;
end

yeo_colormap_glasser=zeros(360,3);
for i=1:360
    key=yeo_in_glasser(i);
    yeo_colormap_glasser(i,:)=yeo_colormap(key,:);

end
figure;
aa=parcel2full((1:361)',glasser_label+1);
aa(isnan(aa))=0;
obj=plot_cortical(aa, 'surface_name','fsa5_inf', 'color_range', [1 361]);
obj.Colormap((colormap([1 1 1;glassercontinouscolors])))

figure
aa=parcel2full((1:7)',yeo_in_glasser);
aa(isnan(aa))=0;
obj=plot_cortical(aa, 'surface_name','fsa5_inf');
obj.Colormap((colormap([1 1 1;yeo_colormap_glasser])))




%% Without combat -

% Retrieve the sample of interest -- either based on group or site
sc_index=strcmp({parameters.group}','SC');
lb_index=strcmp({parameters.group}','LB');
cb_index=strcmp({parameters.group}','CB');
parameters_subset=parameters(logical(sc_index+cb_index));
parameters_subset(contains({parameters_subset.id}','sclb')) = []; % remove sclb of canada
 
% Demographics
age=[parameters_subset.age]';
sex={parameters_subset.sex}';
group={parameters_subset.group}';

mean(age(strcmp(group,'CB')))
std(age(strcmp(group,'CB')))

mean(age(strcmp(group,'SC')))
std(age(strcmp(group,'SC')))
[h,p,cc,stats]=ttest2(age(strcmp(group,'SC')),age(strcmp(group,'CB')))
[tbl,chi2,p] = crosstab(sex,group)

%
gradients_all=zeros(size(parameters_subset,2), size(parameters_subset(1).gradients_sp80_na,1),size(parameters_subset(1).gradients_sp80_na,2));
for i=1:size(parameters_subset,2)
    gradients_all(i,:,:)=parameters_subset(i).gradients_sp90_cs;
end
figure
for i=1:10
    subplot(2,5,i)
    imagesc(corr(gradients_all(:,:,i)'))
    colorbar
    title(i)
end



% Get the covars and residuals
sex=dummyvar(grp2idx({parameters_subset.sex}'));
sex=sex(:,1);
age=[parameters_subset.age]';
site=dummyvar(grp2idx({parameters_subset.site}'));
site=site(:,(1:size(site,2)-1));
covars= [ones(size(gradients_all,1),1) age sex];

residuals=zeros(size(gradients_all));
for i=1:size(gradients_all,2)
    for j=1:size(gradients_all,3)
        [~,~,r] = regress(gradients_all(:,i,j),covars);
        residuals(:,i,j)=r;
    end
end
controlXcb_p=zeros(360,10);controlXcb_t=zeros(360,10);
%Group-level analyses
for i=1:size(gradients_all,2)
    disp(i)
    for j=1:size(gradients_all,3)
        [~,p,~,stats]=ttest2(residuals(strcmp({parameters_subset.group}','CB'),i,j),residuals(strcmp({parameters_subset.group}','SC'),i,j));
        controlXcb_p(i,j)=p;
        controlXcb_t(i,j)=stats.tstat;
        % controlXcb_d(i,j)=abs(computeCohen_d(squeeze(residuals_cbset(i,j,1:14)),squeeze(residuals_cbset(i,j,15:end))));
    end
end
figure
imagesc(fdr_bh(controlXcb_p).*controlXcb_t)
xlabel('Gradient Number')
ylabel('ROI number')
yticks([1:360]')
yticklabels(glasser_metadata.regionName)
% figure
% imagesc(fdr_bh(cbXlb_p).*cbXlb_t)
fdr_mask=(fdr_bh(controlXcb_p).*controlXcb_t);
fdr_mask(isnan(fdr_mask))=0;
obj=plot_hemispheres(fdr_mask(:,1:4), {surf_lh_pial  ,surf_rh_pial}, ...
    'parcellation',glasser_label, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3','Gradient 4'});
obj.colormaps(slanCM('seismic'))


% Decompose significant regions and check the visual areas

glasser_metadata_copy=glasser_metadata;

ix=logical(fdr_mask(:,1:3));
ix=logical(sum(ix,2));
sig_regions=glasser_metadata(ix,:);
sig_regions.g1_ts=fdr_mask((ix),1);
sig_regions.g1_ps=controlXcb_p(ix,1);
sig_regions.g2_ts=fdr_mask((ix),2);
sig_regions.g2_ps=controlXcb_p(ix,2);
sig_regions.g3_ts=fdr_mask((ix),3);
sig_regions.g3_ps=controlXcb_p(ix,3);
sig_regions = sortrows(sig_regions,"g2_ts","descend");

%% Difference in explained variances
lambdas = [parameters_subset.lambdas_sp90_cs].';
mean(lambdas).*100

std(lambdas)
clear controlXcb_t controlXcb_p
residuals=zeros(size(lambdas));
for i=1:size(residuals,2)
    [~,~,r] = regress(lambdas(:,i),covars);
    residuals(:,i)=r;
end
for i=1:size(lambdas,2)
    disp(i)
    [~,p,~,stats]=ttest2(residuals(strcmp({parameters_subset.group}','CB'),i),...
        residuals(strcmp({parameters_subset.group}','SC'),i));
    controlXcb_p(i)=p;
    controlXcb_t(i)=stats.tstat;
    % controlXcb_d(i,j)=abs(computeCohen_d(squeeze(residuals_cbset(i,j,1:14)),squeeze(residuals_cbset(i,j,15:end))));

end

% variances
clear controlXcb_t controlXcb_p
vars=squeeze(var(gradients_all,[],2));
residuals=zeros(size(vars));
for i=1:size(residuals,2)
    [~,~,r] = regress(vars(:,i),covars);
    residuals(:,i)=r;
end
for i=1:size(vars,2)
    disp(i)
    [~,p,~,stats]=ttest2(residuals(strcmp({parameters_subset.group}','CB'),i),...
        residuals(strcmp({parameters_subset.group}','SC'),i));
    controlXcb_p(i)=p;
    controlXcb_t(i)=stats.tstat;
    % controlXcb_d(i,j)=abs(computeCohen_d(squeeze(residuals_cbset(i,j,1:14)),squeeze(residuals_cbset(i,j,15:end))));

end
fdr_bh(controlXcb_p).*controlXcb_t
%% Range

for i=1:size(gradients_all,1)
    for j=1:10
        range(i,j)=squeeze(max(gradients_all(i,:,j))-min(gradients_all(i,:,j)));
    end
end
residuals=zeros(size(range));
for i=1:size(range,2)
    [~,~,r] = regress(range(:,i),[ones(size(range,1),1) age sex]);
    residuals(:,i)=r;
end
clear controlXcb_t controlXcb_p
for i=1:10

    [h,p,cc,stats]=ttest2(residuals(strcmp({parameters_subset.group}','CB'),i),...
        residuals(strcmp({parameters_subset.group}','SC'),i));
    controlXcb_p(i)=p;
    controlXcb_t(i)=stats.tstat;
end
fdr_bh(controlXcb_p)

% Correlation between range and lambdas
[a,b]=corr(lambdas(strcmp({parameters_subset.group}','SC'),1),range(strcmp({parameters_subset.group}','SC'),1),'Type',"Spearman")
[a,b]=corr(lambdas(strcmp({parameters_subset.group}','SC'),2),range(strcmp({parameters_subset.group}','SC'),2),'Type',"Spearman")
[a,b]=corr(lambdas(strcmp({parameters_subset.group}','SC'),3),range(strcmp({parameters_subset.group}','SC'),3),'Type',"Spearman")


%% Rank and distance

% residuals=zeros(size(gradients_all));
% for i=1:size(gradients_all,2)
%     for j=1:size(gradients_all,3)
%         [~,~,r] = regress(gradients_all(:,i,j),covars);
%         residuals(:,i,j)=r;
%     end
% end

% Get ranks
for  i=1:size(parameters_subset,2)
    for j=1:10
        values=[gradients_all(i,strcmp(glasser_metadata.regionName,'V1_L'),j) gradients_all(i,strcmp(glasser_metadata.regionName,'V1_R'),j);...
            gradients_all(i,strcmp(glasser_metadata.regionName,'V2_L'),j) gradients_all(i,strcmp(glasser_metadata.regionName,'V2_R'),j);...
            gradients_all(i,strcmp(glasser_metadata.regionName,'V3_L'),j) gradients_all(i,strcmp(glasser_metadata.regionName,'V3_R'),j);...
            gradients_all(i,strcmp(glasser_metadata.regionName,'V4_L'),j) gradients_all(i,strcmp(glasser_metadata.regionName,'V4_R'),j)];

        values=mean(values,2);
        values=zscale(values,1,0);
        [~,ii]=sort(values,'descend');
        rank = zeros(size(values));
        rank(ii) = 1:length(values);
        ranks(i,:,j)=rank';
    end
end

contingency_tables=cell(size(ranks,2), 10);


for i=1:size(ranks,2)
    for j=1:10

        SC_group=(ranks(strcmp({parameters_subset.group}','SC'),i,j));
        CB_group=(ranks(strcmp({parameters_subset.group}','CB'),i,j));


        [p,h,stats] = ranksum(CB_group,SC_group);
        ranks_p(i,j)=p;
        ranks_z(i,j)=stats.zval;

        modes_sc(i,j) = mode(SC_group);
        modes_cb(i,j) = mode(CB_group);

        medians_sc(i,j) = median(SC_group);
        medians_cb(i,j) = median(CB_group);

        % Get frequencies of values 1, 2, 3, 4 for both groups
        SC_freq = histcounts(SC_group, 0.5:4.5);
        CB_freq = histcounts(CB_group, 0.5:4.5);

        contingency_tables{i,j} = [SC_freq; CB_freq]';
    end
end
fdr_bh(ranks_p(:,1:3)).*ranks_z(:,1:3)
tbl=table(ranks(:,:,1),{parameters_subset.group}');
tbl=sortrows(tbl,"Var2","descend");
imagesc(tbl.Var1')

%================ Min Distance
for  i=1:size(gradients_all,1)
    for j=1:10
        values=[gradients_all(i,strcmp(glasser_metadata.regionName,'V1_L'),j) gradients_all(i,strcmp(glasser_metadata.regionName,'V1_R'),j);...
            gradients_all(i,strcmp(glasser_metadata.regionName,'V2_L'),j) gradients_all(i,strcmp(glasser_metadata.regionName,'V2_R'),j);...
            gradients_all(i,strcmp(glasser_metadata.regionName,'V3_L'),j) gradients_all(i,strcmp(glasser_metadata.regionName,'V3_R'),j);...
            gradients_all(i,strcmp(glasser_metadata.regionName,'V4_L'),j) gradients_all(i,strcmp(glasser_metadata.regionName,'V4_R'),j)];

        values=mean(values,2);
        for k=1:size(values,1)
            difs=values-values(k);
            difs(difs==0)=[];
            distance(k)=(min(abs(difs)));
        end
        min_distances(i,:,j)=distance;
     end
end
% imagesc(min_distances(:,1,1)')
for  i=1:size(min_distances,2)
    for j=1:10
        [h,p,cc,stats]=ttest2(min_distances(strcmp(group,'CB'),i,j),...
            min_distances(strcmp(group,'SC'),i,j));
        distances_t(i,j)=stats.tstat;
        distances_p(i,j)=p;
    end
end
fdr_bh(distances_p(:,1:3)).*distances_t(:,1:3)
%% Structural measures

volume = [parameters_subset.volume]';
thickness = [parameters_subset.thickness_fs].';
surface_area = [parameters_subset.surface_area]';
sulcal_depth = [parameters_subset.sulcal_depth]';
curvature_roi_mp = [parameters_subset.curvature_roi_mp]';



% comparison at group level
covars=[ones(size(thickness,1),1) age sex];
 
for i=1:size(thickness,2)
    [~,~,r] = regress(volume(:,i),covars);
    residuals_vl(:,i)=r;

    [~,~,r] = regress(thickness(:,i),covars);
    residuals_th(:,i)=r;

    [~,~,r] = regress(curvature_roi_mp(:,i),covars);
    residuals_curv(:,i)=r;

    [~,~,r] = regress(surface_area(:,i),covars);
    residuals_sa(:,i)=r;

    [~,~,r] = regress(sulcal_depth(:,i),covars);
    residuals_sd(:,i)=r;
end
clear controlXcb_t controlXcb_p
for i=1:size(thickness,2)

    [h,p,cc,stats]=ttest2(residuals_sd(strcmp({parameters_subset.group}','CB'),i),...
        residuals_sd(strcmp({parameters_subset.group}','SC'),i));
    controlXcb_p_sd(i)=p;
    controlXcb_t_sd(i)=stats.tstat;

    [h,p,cc,stats]=ttest2(residuals_vl(strcmp({parameters_subset.group}','CB'),i),...
        residuals_vl(strcmp({parameters_subset.group}','SC'),i));
    controlXcb_p_vl(i)=p;
    controlXcb_t_vl(i)=stats.tstat;

    [h,p,cc,stats]=ttest2(residuals_th(strcmp({parameters_subset.group}','CB'),i),...
        residuals_th(strcmp({parameters_subset.group}','SC'),i));
    controlXcb_p_th(i)=p;
    controlXcb_t_th(i)=stats.tstat;

    [h,p,cc,stats]=ttest2(residuals_curv(strcmp({parameters_subset.group}','CB'),i),...
        residuals_curv(strcmp({parameters_subset.group}','SC'),i));
    controlXcb_p_ci(i)=p;
    controlXcb_t_ci(i)=stats.tstat;

    [h,p,cc,stats]=ttest2(residuals_sa(strcmp({parameters_subset.group}','CB'),i),...
        residuals_sa(strcmp({parameters_subset.group}','SC'),i));
    controlXcb_p_sa(i)=p;
    controlXcb_t_sa(i)=stats.tstat;
end
plot_hemispheres([(fdr_bh(controlXcb_p_vl).*controlXcb_t_vl)' ...
    (fdr_bh(controlXcb_p_th).*controlXcb_t_th)'...
    (fdr_bh(controlXcb_p_sa).*controlXcb_t_sa)'...
    (fdr_bh(controlXcb_p_sd).*controlXcb_t_sd)'...
    (fdr_bh(controlXcb_p_ci).*controlXcb_t_ci)'],...
    {surf_lh,surf_rh},'parcellation',glasser_label,'labeltext',{ 'Volume', 'Thickness','Surface Area','Sulcal Depth', 'Curvature'});


for i=1:size(parameters_subset,2)

    [coeff,score,latent,tsquared,explained] = pca([volume(i,:)' thickness(i,:)' surface_area(i,:)'  curvature_roi_mp(i,:)' sulcal_depth(i,:)']);
    %
    pca_str(i,:)=score(:,1);
    pca_var(i,:)=explained;
end

plot_hemispheres([mean(pca_str(strcmp({parameters_subset.group}','SC'),:))'...
    mean(pca_str(strcmp({parameters_subset.group}','CB'),:))'], {surf_lh,surf_rh},'parcellation',glasser_label);


clear controlXcb_t controlXcb_p
for i=1:size(pca_str,2)

    [h,p,cc,stats]=ttest2(pca_str(strcmp({parameters_subset.group}','CB'),i),...
        pca_str(strcmp({parameters_subset.group}','SC'),i));
    controlXcb_p(i)=p;
    controlXcb_t(i)=stats.tstat;
end
plot_hemispheres(((controlXcb_p<0.05).*controlXcb_t)', {surf_lh,surf_rh},'parcellation',glasser_label);
%
str_results=[((controlXcb_p<0.05).*controlXcb_t)' (fdr_bh(controlXcb_p_vl).*controlXcb_t_vl)' ...
    (fdr_bh(controlXcb_p_th).*controlXcb_t_th)'...
    (fdr_bh(controlXcb_p_sa).*controlXcb_t_sa)'...
    (fdr_bh(controlXcb_p_sd).*controlXcb_t_sd)'...
    (fdr_bh(controlXcb_p_ci).*controlXcb_t_ci)']


residuals=zeros(size(gradients_all(1:71,:,:)));
covars=[ones(size(residuals,1),1) age(1:71) sex(1:71)];
for i=1:size(residuals,2)
    for j=1:10
        [~,~,r] = regress(gradients_all(1:71,i,j),covars);
        residuals(:,i,j)=r;
    end
end
parameters_2=parameters_subset(1:71);

batch=(grp2idx({parameters_2.site}'));
structural_harmonized = combat(pca_str(1:71,:)', batch, [age(1:71) sex(1:71)], 1)';

for i=1:360
    for j=1:10
        xx(i,j)=corr(residuals(strcmp({parameters_2.group}','CB'),i,j),structural_harmonized(strcmp({parameters_2.group}','CB'),j));
        yy(i,j)=corr(residuals(strcmp({parameters_2.group}','SC'),i,j),structural_harmonized(strcmp({parameters_2.group}','SC'),j));
    end
end
for i=1:360
    for j=1:3
        r1 = xx(i,j); % Correlation coefficient for group 1
        r2 = yy(i,j); % Correlation coefficient for group 2
        n1 = 27; % Sample size for group 1
        n2 = 44;
        % Fisher Z-transformation
        z1 = 0.5 * log((1 + r1) / (1 - r1));
        z2 = 0.5 * log((1 + r2) / (1 - r2));
        % Calculate the Z difference
        z_diff = (z1 - z2) / sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)));
        % P-value for two-tailed test
        pvals(i,j)=2 * (1 - normcdf(abs(z_diff)));
    end
end
plot_hemispheres([xx(:,1) yy(:,1) (xx(:,1)-yy(:,1))], {surf_lh,surf_rh},'parcellation',glasser_label,'labeltext', {'CB','SC', 'Difference'});
plot_hemispheres([xx(:,2) yy(:,2) (xx(:,2)-yy(:,2))], {surf_lh,surf_rh},'parcellation',glasser_label,'labeltext', {'CB','SC', 'Difference'});
plot_hemispheres([xx(:,3) yy(:,3) (xx(:,3)-yy(:,3))], {surf_lh,surf_rh},'parcellation',glasser_label,'labeltext', {'CB','SC', 'Difference'});



mean(abs(xx(:,1:3)))
std(abs(xx(:,1:3)))

mean(abs(yy(:,1:3)))
std(abs(yy(:,1:3)))


[h,p,~,stats]=ttest2(xx(:,1),xx(:,2))

[h,p,~,stats]=ttest2(xx(:,1),yy(:,1))

[h,p,~,stats]=ttest2(yy(:,1),yy(:,2))

[h,p,~,stats]=ttest2(xx(:,2),yy(:,2))

[h,p,~,stats]=ttest2(xx(:,3),yy(:,3))
[h,p,~,stats]=ttest2(yy(:,2),yy(:,3))

[a,b]=corr(xx(:,1),mean(gradients_all(strcmp({parameters_subset.group}','SC'),:,1))','type','Spearman')
[a,b]=corr(yy(:,1),mean(gradients_all(strcmp({parameters_subset.group}','CB'),:,1))','type','Spearman')
[a,b]=corr(xx(:,2),mean(gradients_all(strcmp({parameters_subset.group}','SC'),:,2))','type','Spearman')
[a,b]=corr(yy(:,2),mean(gradients_all(strcmp({parameters_subset.group}','CB'),:,2))','type','Spearman')
[a,b]=corr(xx(:,3),mean(gradients_all(strcmp({parameters_subset.group}','SC'),:,3))','type','Spearman')
[a,b]=corr(yy(:,3),mean(gradients_all(strcmp({parameters_subset.group}','CB'),:,2))','type','Spearman')
% mean_structural_sc=mean(structural_harmonized(strcmp({parameters_2.group}','SC'),:))';
