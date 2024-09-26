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
% cors=corr(gradients_all(:,:,1)');
%
% mean(mean(cors(1:47,1:47)))
% mean(mean(cors(48:107,48:107)))
% mean(mean(cors(107:end,107:end)))
% [h,p,cc,stats]=ttest2(mean(cors(1:47,1:47)),mean(cors(1:47,48:107)))
%
% mean(mean(cors(1:47,48:107)))
% mean(mean(cors(1:47,107:end)))
% mean(mean(cors(48:107,48:107)))


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
%% Load the full time series


dset='krakow';

ts_clean=load(['/media/koba/MULTIBOOT/blindness_gradients/datasets/' dset '/derivatives/ts_clean_' dset '.mat']);
ts_clean=ts_clean.ts_clean;
participants=readtable(['/media/koba/MULTIBOOT/blindness_gradients/datasets/' dset '/participants.csv']);


ids={parameters_subset(strcmp({parameters_subset.site}',dset)).id}';
ids=erase(ids,[dset '_']);
mask=zeros(size(participants,1),1);
for i=1:size(ids,1)
    mask=mask+contains(participants.participant_id,ids{i});
end
ts_clean=ts_clean(logical(mask),:,:);

ts_clean3=ts_clean;
% ts_clean_vertexwise=[ts_clean1;ts_clean2;ts_clean3];
ts_clean_network=ts_clean1(:,logical(network_regions_vertex),:);
network_matrices1=zeros(size(ts_clean_network,1),size(ts_clean_network,2),size(ts_clean_network,2));
for i=1:size(ts_clean_network,1)
    disp(i)
    network_matrix=atanh(corr(squeeze(ts_clean_network(i,:,:))'));
    network_matrix(isinf(network_matrix))=1;
    network_matrices1(i,:,:)=network_matrix;
end
ts_clean_network=ts_clean2(:,logical(network_regions_vertex),:);
network_matrices2=zeros(size(ts_clean_network,1),size(ts_clean_network,2),size(ts_clean_network,2));
for i=1:size(ts_clean_network,1)
    disp(i)
    network_matrix=atanh(corr(squeeze(ts_clean_network(i,:,:))'));
    network_matrix(isinf(network_matrix))=1;
    network_matrices2(i,:,:)=network_matrix;
end
ts_clean_network=ts_clean3(:,logical(network_regions_vertex),:);
network_matrices3=zeros(size(ts_clean_network,1),size(ts_clean_network,2),size(ts_clean_network,2));
for i=1:size(ts_clean_network,1)
    disp(i)
    network_matrix=atanh(corr(squeeze(ts_clean_network(i,:,:))'));
    network_matrix(isinf(network_matrix))=1;
    network_matrices3(i,:,:)=network_matrix;
end
network_matrices=[network_matrices1;network_matrices2;network_matrices3];

clear  network_matrices1 network_matrices2 network_matrices3 ts_clean_network network_matrix gradients_all

%%
reference_network= GradientMaps('kernel','cs','approach','dm','alignment','pa','n_components',10);
mean_reference=squeeze(mean(network_matrices(strcmp({parameters_subset.group}','SC'),:,:)));
reference_network = reference_network.fit(mean_reference, 'sparsity',90);
mean_reference_network_full=zeros(size(reference_gradients,1),size(reference_network.aligned{1},2));
for i=1:size(reference_network.aligned{1},2)
    mean_reference_network_full(network_regions_vertex==1,i)=reference_network.aligned{1}(:,i);
end
plot_hemispheres(mean_reference_network_full(:,1:3), {surf_lh,surf_rh});
gradient_in_euclidean(mean_reference_network_full(:,1:2),{surf_lh,surf_rh});
scree_plot(reference_network.lambda{1})

for i=1:size(network_regions,1)
    ix=glasser_label==network_regions(i);
    means_network(i,:)=mean(mean_reference_network_full(ix,:));
    sds_network(i,:)=std(mean_reference_network_full(ix,:));
end
network_colors=grp2idx(glasser_metadata.cortex(network_regions,:))/max(grp2idx(glasser_metadata.cortex(network_regions,:)));
% network_colors=[network_colors network_colors/2 network_colors*2];


for j=1:3
    subplot(3,1,j)
    tbl=table(means_network(:,j),glasser_metadata(glasser_metadata.network_index==1,:).regionName,network_colors);
    tbl=sortrows(tbl,'Var1','descend');
    b=bar(tbl.Var1,'FaceColor','flat','CData',tbl.Var3);
    xticks([1:size(network_regions,1)])
    xticklabels(tbl.Var2)
    colororder("earth")
    xlabel('Regions')
    ylabel('Gradient scores')
    title({'Gradient ' string(j)})
end


gradients_network_all=zeros(size(network_matrices,1),size(network_matrices,2),size(reference_gradients,2));
lambdas=zeros(size(reference_network.lambda{1},1),size(reference_gradients,2));
for i=1:size(network_matrices,1)
    disp(i)
        m=(squeeze(network_matrices(i,:,:)));
        gradfit=GradientMaps('kernel','cs','approach','dm','alignment','pa');
        gradfit=gradfit.fit(m,'reference',reference_network.gradients{1},'sparsity',90);
        gradients_network_all(i,:,:)=gradfit.aligned{1};
        lambdas(:,i)=gradfit.lambda{1}./sum(gradfit.lambda{1});
    
end


% Continue here after loading the gradients 
sex=dummyvar(grp2idx({parameters_subset.sex}'));
sex=sex(:,1);
age=[parameters_subset.age]';
site=dummyvar(grp2idx({parameters_subset.site}'));
site=site(:,(1:size(site,2)-1));
covars= [ones(size(gradients_network_all,1),1) age sex];

residuals=zeros(size(gradients_network_all));
for i=1:size(gradients_network_all,2)
    disp(i)
    for j=1:size(gradients_network_all,3)
        [~,~,r] = regress(gradients_network_all(:,i,j),covars);
        residuals(:,i,j)=r;
    end
end


controlXcb_p=zeros(360,10);controlXcb_t=zeros(360,10);
%Group-level analyses
for i=1:size(gradients_network_all,2)
    disp(i)
    for j=1:size(gradients_network_all,3)
        [~,p,~,stats]=ttest2(residuals(strcmp({parameters_subset.group}','CB'),i,j),residuals(strcmp({parameters_subset.group}','SC'),i,j));
        controlXcb_p(i,j)=p;
        controlXcb_t(i,j)=stats.tstat;
        % controlXcb_d(i,j)=abs(computeCohen_d(squeeze(residuals_cbset(i,j,1:14)),squeeze(residuals_cbset(i,j,15:end))));
    end
end
figure
imagesc(fdr_bh(controlXcb_p).*controlXcb_t)
xlabel('Gradient number ')
ylabel('Vertex')
title('Significant differences in gradient scores between SC and CB')
fdr_mask=(fdr_bh(controlXcb_p).*controlXcb_t);

mean_network_full_diff=zeros(size(network_regions_vertex,1),size(gradients_network_all,3));

for j=1:10
    mean_network_full_diff(network_regions_vertex==1,j)=fdr_mask(:,j);
end


% figure
% imagesc(fdr_bh(cbXlb_p).*cbXlb_t)
obj=plot_hemispheres(mean_network_full_diff(:,1:4), {surf_lh  ,surf_rh}, ...
    'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3','Gradient 4'});  


mean_sc=squeeze(mean(gradients_network_all(strcmp({parameters_subset.group}','SC'),:,:)));
mean_cb=squeeze(mean(gradients_network_all(strcmp({parameters_subset.group}','CB'),:,:)));
%mean=[mean_sc(:,1:2) mean_cb(:,1:2)];
mean_full_sc=zeros(size(network_regions_vertex,1),size(gradients_network_all,3));
mean_full_cb=zeros(size(network_regions_vertex,1),size(gradients_network_all,3));


for j=1:4
    mean_full_cb(network_regions_vertex==1,j)=mean_cb(:,j);
    mean_full_sc(network_regions_vertex==1,j)=mean_sc(:,j);
end
plot_hemispheres([mean_full_sc(:,1:2) mean_full_cb(:,1:2)], {surf_lh  ,surf_rh});  

gradient_in_euclidean(mean_sc(:,1:2),{surf_lh,surf_rh});
gradient_in_euclidean(mean_cb(:,1:2),{surf_lh,surf_rh});

residuals=zeros(size(lambdas));
for i=1:size(residuals,2)
    [~,~,r] = regress(lambdas(:,i),covars);
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



for i=1:size(gradients_all,1)
    for j=1:10
        range(i,j)=squeeze(max(gradients_network_all(i,:,j))-min(gradients_network_all(i,:,j)));
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


clear controlXcb_t controlXcb_p
vars=squeeze(var(gradients_network_all,[],2));
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
% Mean over networks


network_regions=find(glasser_metadata.network_index==1);
unique_region_names=unique(glasser_metadata.region(network_regions),'stable');
zzz=glasser_metadata.cortex(network_regions);
cortex_names=zzz(1:26);
for i=1:size(unique_region_names,1)
    region_codes=find(strcmp(glasser_metadata.region,unique_region_names(i)));
    ix=ismember(glasser_label,region_codes);
    means_network(i,:)=mean(to_draw(ix,:));
    sds_network(i,:)=std(to_draw(ix,:));
end


fff=mean_full_sc(:,1);

%%
thickness = [parameters_subset.thickness_mp].';
curvature = [parameters_subset.curvature_mp].';
sulcal_depth = [parameters_subset.sulcal_depth_vertex].';
surface_area = [parameters_subset.area_vertex].';
volume = [parameters_subset.volume_vertex].';

for i=1:size(parameters_subset,2)
    % structural_matrices(i,:,:) = [volume(i,:)' thickness(i,:)' curvature_int(i,:)' curvature_rect(i,:)' curvature_rect_gauss(i,:)' surface_area(i,:)' folding_index(i,:)'];
    % structural_conn_matrices(i,:,:)=corr(squeeze(structural_matrices(i,:,:))');
    [coeff,score,latent,tsquared,explained] = pca([volume(i,:)' thickness(i,:)' curvature(i,:)' surface_area(i,:)' sulcal_depth(i,:)']);
    %
    pca_str(i,:)=score(:,1);
    pca_var(i,:)=explained;
end
mean(pca_var)
std(pca_var)
plot_hemispheres([mean(pca_str(strcmp({parameters_subset.group}','SC'),:))'...
    mean(pca_str(strcmp({parameters_subset.group}','CB'),:))'], {surf_lh,surf_rh});

str_vertex=[mean(pca_str(strcmp({parameters_subset.group}','SC'),:))'...
    mean(pca_str(strcmp({parameters_subset.group}','CB'),:))'];
% for i=1:size(thickness,2)
%     [~,~,r] = regress(pca_str(:,i),covars);
%     residuals(:,i)=r;
% end
% 
% for i=1:size(thickness,2)
% 
%     [h,p,cc,stats]=ttest2(residuals(strcmp({parameters_subset.group}','CB'),i),...
%         residuals(strcmp({parameters_subset.group}','SC'),i));
%     controlXcb_p_th(i)=p;
%     controlXcb_t_th(i)=stats.tstat;
% 
% end
% 
% plot_hemispheres((fdr_bh(controlXcb_p_th).*controlXcb_t_th)',{surf_lh,surf_rh})


parameters_2=parameters_subset(1:71);
batch=(grp2idx({parameters_2.site}'));
% nonzero_ix=sum(thickness~=0)==85;
% thicknes_nonzero=thickness(:,nonzero_ix);
structural_network=(pca_str(:,logical(network_regions_vertex)));


structural_harmonized = combat(structural_network(1:71,:)', batch, [age(1:71) sex(1:71)], 1)';
residuals=zeros(size(gradients_network_all(1:71,:,:)));
for i=1:size(gradients_network_all,2)
    for j=1:3
        [~,~,r] = regress(gradients_network_all(1:71,i,j),covars(1:71,:));
        residuals(:,i,j)=r;
    end
end
%% Try an approach that will smooth the data
glasser_metadata.network_index=(strcmp(glasser_metadata.cortex,'Primary_Visual') +...
    strcmp(glasser_metadata.cortex,'Early_Visual') + ...
    strcmp(glasser_metadata.cortex,'Dorsal_Stream_Visual') +...
    strcmp(glasser_metadata.cortex,'Ventral_Stream_Visual')+...
    strcmp(glasser_metadata.cortex,'MT+_Complex_and_Neighboring_Visual_Areas') );


network_regions=find(glasser_metadata.network_index==1);
network_regions_vertex=zeros(size(glasser_label,1),1);
for i=1:length(network_regions)
    network_regions_vertex=double(glasser_label==network_regions(i))+network_regions_vertex;
end
    plot_hemispheres(double(network_regions_vertex), {surf_lh,surf_rh}   );

% Get the smaller regions in the network
parameters_2=parameters_subset(1:71);
regions_of_interest={'V1','V2','V3','V4'};
for k=1:length(regions_of_interest)
    glasser_metadata.single_region_index=(strcmp(glasser_metadata.region,regions_of_interest{k}));
    single_region=find(glasser_metadata.single_region_index==1);
    single_region_vertex=zeros(size(glasser_label,1),1);
    for i=1:length(single_region)
        single_region_vertex=double(glasser_label==single_region(i))+single_region_vertex;
    end
    % plot_hemispheres(double(single_region_vertex), {surf_lh,surf_rh}   );

    single_region_in_network_ix=logical(single_region_vertex(logical(network_regions_vertex)));

    for i=1:size(parameters_2,2)
        for j=1:3
            [a,b]=corr(residuals(i,single_region_in_network_ix,j)',structural_harmonized(i,single_region_in_network_ix)','Type','Spearman');
            corrs(k,i,j)=a;
            corrs_p(k,i,j)=b;
        end
    end
end

sc_rgb=hex2rgb('#440154');
cb_rgb=hex2rgb('#43BF71');

colors = [sc_rgb';cb_rgb'];

cnt = 0;
for i = 1:size(corrs, 1)
    for j = 1:size(corrs, 3)
        cnt = cnt + 1;
        subplot(4, 3, cnt)
        
        boxC = boxchart(corrs(i, :, j), 'GroupByColor', {parameters_2.group}');
        
       boxC(1).BoxFaceColor = colors(2, :);  
        boxC(2).BoxFaceColor = colors(1, :);
        if i==3
        legend
        end
        xlabel(regions_of_interest(i))
        ylabel('Rho')
    end
end
boxC = boxchart(corrs(1, :, 2), 'GroupByColor', {parameters_2.group}');
        
       boxC(1).BoxFaceColor = colors(2, :);  
        boxC(2).BoxFaceColor = colors(1, :);
        if i==3
        legend
        end
        xlabel(regions_of_interest(i))
        ylabel(['Gradient ' + string(j)])

mean(corrs(4,strcmp({parameters_2.group}','SC'),2))
std(corrs(4,strcmp({parameters_2.group}','SC'),2))

mean(corrs(4,strcmp({parameters_2.group}','CB'),2))
std(corrs(4,strcmp({parameters_2.group}','CB'),2))
