clear
parameters1=load('/media/koba/MULTIBOOT/blindness_gradients/datasets/canada/derivatives/parameters.mat');
parameters1=parameters1.parameters;
parameters2=load('/media/koba/MULTIBOOT/blindness_gradients/datasets/baltimore/derivatives/parameters.mat');
parameters2=parameters2.parameters;
parameters3=load('/media/koba/MULTIBOOT/blindness_gradients/datasets/krakow/derivatives/parameters.mat');
parameters3=parameters3.parameters;

parameters=[parameters1, parameters2, parameters3];
% parameters(111).gradients=zeros(360,10); parameters(111).lambdas=zeros(17,1);

% glasser_label=table2array(readtable('/media/koba/MULTIBOOT/blindness_gradients/source/parcellations/glasser_fsaverage5_labels.txt'));
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

glasser_metadata=readtable('/media/koba/MULTIBOOT/blindness_gradients/source/parcellations/HCP-MMP1_UniqueRegionList.csv');
glasser5networks=table2array(readtable('/media/koba/MULTIBOOT/blindness_gradients/source/parcellations/glasser5networks.txt'));
glasser5networks=glasser5networks(1:360);
glassercontinouscolors=table2array(readtable('/media/koba/MULTIBOOT/blindness_gradients/source/parcellations/glasser_continous_colors.txt'));
glassercontinouscolors=glassercontinouscolors(:,1:3);
label=glasser_label;
addpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/BrainStat-master'))
reference_gradients=fetch_gradients();
[surf_lh, surf_rh] = fetch_template_surface('fsaverage5','layer','inflated');
rmpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/BrainStat-master'))



age = ([parameters.age].');
sex = ([parameters.sex].');
group = {parameters.group}.';
parenchyma = ([parameters.parenchyma].');
gm = [parameters.gm].';
wm = [parameters.wm].';
subcorticals=[parameters.L_Thal; parameters.L_Caud; parameters.L_Puta; parameters.L_Pall; parameters.L_Hippo; parameters.L_Amyg; parameters.L_Acc; parameters.R_Thal; parameters.R_Caud; parameters.R_Puta; parameters.R_Pall; parameters.R_Hippo; parameters.R_Amyg; parameters.R_Acc; parameters.Brain_Stem].';;
thickness_mp = [parameters.thickness_mp].';
thickness_roi_mp = [parameters.thickness_roi_mp].';
curvature_mp = [parameters.curvature_mp].';
curvature_roi_mp = [parameters.curvature_roi_mp].';
thickness_fs = [parameters.thickness_fs].';
folding_index = [parameters.folding_index].';
rectified_curvature = [parameters.rectified_curvature].';
rectified_curvature_gaussian = [parameters.rectified_curvature_gaussian].';
intrinsic_curvature = [parameters.intrinsic_curvature].';




lambdas = [parameters.lambdas_sp90_cs].';

for i=1:size(parameters,2)
    corrmats(i,:,:)=parameters(i).corrmat;
end
corrmats_mean=mean(corrmats,3);

gd=[parameters.gd];
for i=1:size(parameters,2)
    gds(i,:,:)=parameters(i).gd;
end
gds_mean=mean(gds,3);

% gradients=[parameters.gradients];
for i=1:size(parameters,2)
    gradients_sp90_cs(i,:,:)=parameters(i).gradients_sp90_cs;
    gradients_sp90_na(i,:,:)=parameters(i).gradients_sp90_na;
    gradients_sp80_na(i,:,:)=parameters(i).gradients_sp80_na;
end

tissue_names={'Parenchymal volume' 'GM volume' 'WM volume' 'L Thal' 'L Caud' 'L Puta' 'L Pall' 'L Hippo' 'L Amyg' 'L Acc' ...
    'R Thal' 'R Caud' 'R Puta' 'R Pall' 'R Hippo' 'R Amyg' 'R Acc' ...
    'Brain Stem'};
tissues=[parenchyma gm wm subcorticals];

for i=1:size(tissues,2)
    subplot(size(tissues,2)/3,3,i)
    p=stem((tissues(:,i)),'filled');
    xlabel(tissue_names{i})
    ylabel('Z score')
     % ylim([-3,3])
end
ix=strcmp({parameters.group}','CB');
for i=1:size(tissues,2)
    subplot(size(tissues,2)/3,3,i)
    boxchart((tissues(ix,i)),'GroupByColor',{parameters(ix).site})
    xlabel(tissue_names{i})
    ylabel('Z score')
    % ylim([-3,3])
end

sex_d=dummyvar(grp2idx(sex));
mod = [age sex_d(:,2)];
data_harmonized = combat(tissues(ix,:)', grp2idx({parameters(ix).site}'), mod(ix,:), 1);

figure
for i=1:size(tissues,2)
    subplot(size(tissues,2)/3,3,i)
    boxchart((data_harmonized(i,:)'),'GroupByColor',{parameters(ix).site})
    xlabel(tissue_names{i})
    ylabel('Z score')
    % ylim([-3,3])
end


subplot(2,4,1)
imagesc(corr(corrmats_mean'))
colorbar
title('Mean Group FC')
subplot(2,4,2)
imagesc(corr(thickness_mp'))
colorbar
title('Mean Thickness')
subplot(2,4,3)
imagesc(corr(thickness_roi_mp'))
colorbar
title('Mean Thickness in Glasser')
subplot(2,4,4)
imagesc(corr(curvature_mp'))
colorbar
title('Mean Curvature')
subplot(2,4,5)
imagesc(corr(curvature_roi_mp'))
colorbar
title('Mean Curvature in Glasser')
subplot(2,4,6)
imagesc(corr(gds_mean'))
colorbar
title('Mean GD')
subplot(2,4,7)
imagesc(corr(rectified_curvature'))
colorbar
title('Rect. Curv.')
subplot(2,4,8)
imagesc(corr(rectified_curvature_gaussian'))
colorbar
title('Rect. Curv. Gaus.')


figure
for i=1:10
    subplot(2,5,i)
    imagesc(corr(gradients_sp90_cs(:,:,i)'))
    colorbar
    title(i)
end



 imagesc(corr(gradients_all(:,:,1)'))
plot_hemispheres(mean(gradients_sp90_cs(:,:,1),'omitmissing')', {surf_lh,surf_rh}, ...
    'parcellation',label, 'labeltext',{'Mean FC'});
gradient_in_euclidean(parameter)
plot_hemispheres(mean(thickness(:,:,1),'omitmissing')', {surf_lh,surf_rh}, ...
    'labeltext',{'Mean FC'});


gradient_in_euclidean(parameters(97).gradients(:,1:2),{surf_lh,surf_rh},glasser_label);

for i=1:size(tissues,2)
    subplot(6,3,i)
    p=stem((lambdas(:,i)),'filled');
    xlabel(tissue_names{i})
    ylabel('Z scored lambda')
    % ylim([-3,3])
end

figure
for i=1:size(tissues,2)
    subplot(6,3,i)
boxchart((lambdas(:,i)),'GroupByColor',{parameters.site})
    xlabel(tissue_names{i})
    ylabel('Lambda')
      % ylim([-3,3])
end

for i=1:size(parameters,2)
    for j=1:10
    gradient_range(i,j)=max(squeeze(gradients_all(i,:,j)))-min(squeeze(gradients_all(i,:,j)));
end
end
figure
for i=1:10
    subplot(5,2,i)
    boxchart(({parameters.site}'),(gradient_range(:,i)),'GroupByColor',{parameters.group})
    % xlabel(tissue_names{i})
    % ylabel('Lambda')
    % ylim([-3,3])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
aaa=[zscore(tissues(1:50,8));
    zscore(tissues(51:107,8));
    zscore(tissues(108:end,8))];

anova1(aaa(1:14),{parameters(51:107).group}')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=age(strcmp({parameters.site}','baltimore'))
b=sex(strcmp({parameters.site}','baltimore'))
c=group(strcmp({parameters.site}','baltimore'))

[h,p,cc,stats]=ttest2(a(strcmp(c,'SC')),a(strcmp(c,'CB')))
[tbl,chi2,p,labels] = crosstab(b,c)

sc_ix=strcmp({parameters.group}','SC')
cb_ix=strcmp({parameters.group}','CB')

[h,p,cc,stats]=ttest2([parameters(sc_ix).age]',[parameters(cb_ix).age]')

[tbl,chi2,p,labels] = crosstab({parameters(logical(sc_ix+cb_ix)).sex}',{parameters(logical(sc_ix+cb_ix)).group}')

