function [ts_clean, parameters]=clean_and_prepare_data(dset_name,TR)

% TR=1.5;dset_name='krakow';
%
% Files related to parcellation
%glasser_label=table2array(readtable('/media/koba/MULTIBOOT/blindness_gradients/source/parcellations/glasser_fsaverage5_labels.txt'));
[~, label_l, colortable_l]=read_annotation('/home/koba/micapipe/parcellations/lh.glasser-360_mics.annot');
[~, label_r, colortable_r]=read_annotation('/home/koba/micapipe/parcellations/rh.glasser-360_mics.annot');
for i=1:size(colortable_l.table,1)
    label_l(label_l==colortable_l.table(i,5))=i;
    label_r(label_r==colortable_r.table(i,5))=i;
end
label_r=label_r+180;
label_r(label_r==181)=1;
glasser_label=[label_l;label_r]-1;

subcortical_labels=[10 11 12 13 17 18 26 49 50 51 52 53 54 58 16];
subcortical_names={'L_Thal' 'L_Caud' 'L_Puta' 'L_Pall' 'L_Hippo' 'L_Amyg' 'L_Acc' ...
    'R_Thal' 'R_Caud' 'R_Puta' 'R_Pall' 'R_Hippo' 'R_Amyg' 'R_Acc' ...
    'Brain_Stem'};


% Parameters for band-pass filtering
lowFreq = 0.01; % Lower cutoff frequency in Hz
highFreq = 0.1; % Upper cutoff frequency in Hz
filterOrder = 4; % Filter order
samplingRate = 1 / TR;
[b, a] = butter(filterOrder, [lowFreq, highFreq] * 2 / samplingRate);


dset_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/';
dset_path=strrep(dset_path,'DSET',dset_name);
fileList = dir(fullfile(dset_path, '**', '*'));
fileList([fileList.isdir].' == 1) = [];



participants_string='participants.csv';
timeseries_string='_hemi-HEMI_surf-fsaverage5.func.gii';
fd_string='_space-func_desc-se_pve_WM.txt';
wm_string='_space-func_desc-se_pve_WM.txt';
csf_string='_space-func_desc-se_pve_CSF.txt';
motion_string='_space-func_desc-se.1D';
thickness_string='_hemi-HEMI_surf-fsaverage5_label-thickness.func.gii';
curvature_string='_hemi-HEMI_surf-fsaverage5_label-curv.func.gii';
gd_string='_atlas-glasser-360_GD.shape.gii';
t1_string='_space-nativepro_T1w_brain_mask.nii.gz';
gm_nii_string='_space-nativepro_T1w_brain_pve_1.nii.gz';
wm_nii_string='_space-nativepro_T1w_brain_pve_2.nii.gz';
first_string='_space-nativepro_T1w_atlas-subcortical.nii.gz';
morphology_string='HEMI.glasser.stats';




participants=readtable([fileList(contains({fileList.name}',participants_string)).folder '/'  fileList(contains({fileList.name}',participants_string)).name]);


%participants_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/participants.csv';
%timeseries_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/func/desc-se_task-rest_bold/surf/SUBID_hemi-HEMI_surf-fsnative.func.gii';
% fd_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/func/desc-se_task-rest_bold/volumetric/SUBID_space-func_desc-se_metric_FD.1D';
% wm_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/func/desc-se_task-rest_bold/volumetric/SUBID_space-func_desc-se_pve_WM.txt';
% csf_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/func/desc-se_task-rest_bold/volumetric/SUBID_space-func_desc-se_pve_CSF.txt';
% motion_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/func/desc-se_task-rest_bold/volumetric/SUBID_space-func_desc-se.1D';
% thickness_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/maps/SUBID_hemi-HEMI_surf-fsaverage5_label-thickness.func.gii';
% curvature_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/maps/SUBID_hemi-HEMI_surf-fsaverage5_label-curv.func.gii';
% gd_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/dist/SUBID_atlas-glasser-360_GD.shape.gii';
% t1_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/anat/SUBID_space-nativepro_T1w_brain_mask.nii.gz';
% gm_nii_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/anat/SUBID_space-nativepro_T1w_brain_pve_1.nii.gz';
% wm_nii_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/anat/SUBID_space-nativepro_T1w_brain_pve_2.nii.gz';
% first_path='/media/koba/MULTIBOOT/blindness_gradients/datasets/DSET/derivatives/micapipe_v0.2.0/SUBID/parc/SUBID_space-nativepro_T1w_atlas-subcortical.nii.gz';
%
% timeseries_path=strrep(timeseries_path,'DSET',dset_name);
% fd_path=strrep(fd_path,'DSET',dset_name);
% wm_path=strrep(wm_path,'DSET',dset_name);
% csf_path=strrep(csf_path,'DSET',dset_name);
% motion_path=strrep(motion_path,'DSET',dset_name);
% thickness_path=strrep(thickness_path,'DSET',dset_name);
% curvature_path=strrep(curvature_path,'DSET',dset_name);
% gd_path=strrep(gd_path,'DSET',dset_name);
% t1_path=strrep(t1_path,'DSET',dset_name);
% gm_nii_path=strrep(gm_nii_path,'DSET',dset_name);
% wm_nii_path=strrep(wm_nii_path,'DSET',dset_name);
% first_path=strrep(first_path,'DSET',dset_name);


% Calculate the reference gradients
addpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/BrainStat-master'))
reference_gradients=fetch_gradients();
rmpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/BrainStat-master'))

for i=1:max(glasser_label)
    for j=1:size(reference_gradients,2)
        reference_gradients_reduced(i,j,:)=mean(squeeze(reference_gradients(glasser_label==i,j,:)),'omitmissing');
    end
end

% Clean the time series
for i=1:size(participants,1)

    disp(['Cleaning the time series of subject ' num2str(i) ])

    timeseries_path_right=strrep(timeseries_string,'HEMI','R');
    lis=[fileList(contains({fileList.name}',timeseries_path_right))];
    timeseries_right=gifti([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    timeseries_right=timeseries_right.cdata;

    timeseries_path_left=strrep(timeseries_string,'HEMI','L');
    lis=[fileList(contains({fileList.name}',timeseries_path_left))];
    timeseries_left=gifti([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    timeseries_left=timeseries_left.cdata;

    timeseries=[timeseries_left;timeseries_right];
    clear timeseries_left timeseries_right

    % Postprocess
    lis=[fileList(contains({fileList.name}',fd_string))];
    fd=table2array(readtable([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name],"FileType","text",'NumHeaderLines',0));

    lis=[fileList(contains({fileList.name}',wm_string))];
    wm=table2array(readtable([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name],"FileType","text",'NumHeaderLines',0));


    lis=[fileList(contains({fileList.name}',csf_string))];
    csf=table2array(readtable([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name],"FileType","text",'NumHeaderLines',0));


    lis=[fileList(contains({fileList.name}',motion_string))];
    motion=table2array(readtable([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name],"FileType","text",'NumHeaderLines',0));

    ortvec=[wm csf motion];

    ortvec_diff=[zeros(1,size(ortvec,2)); diff(ortvec)];
    ortvec_power=ortvec.^2;
    ortvec_diff_power=ortvec_diff.^2;
    linear_trend=[1:size(timeseries,2)]';
    confounds=[fd ortvec ortvec_diff ortvec_power ortvec_diff_power linear_trend];

    timeseries_clean=zeros(size(timeseries));
    for j=1:size(timeseries,1)

        [~,~,residuals] = regress(timeseries(j,:)',[ones(size(timeseries,2),1) confounds]);
        timeseries_clean(j,:) = filter(b, a, residuals);
    end
    ts_clean(i,:,:)=timeseries_clean;

end

% Parcellate the time series
ts_clean_reduced=zeros(size(ts_clean,1),max(glasser_label),size(ts_clean,3));
for i=1:size(participants,1)
    disp(['Reducing the time series for subject ' num2str(i)])
    for j=1:max(glasser_label)
        ts_clean_reduced(i,j,:)=mean(squeeze(ts_clean(i,glasser_label==j,:)),'omitnan');
    end
end



% Extract the relevant measures
for i=1:size(participants,1)
    disp(['Calculating the parameters for subject ' num2str(i)])
    % Demographics
    parameters(i).id=[dset_name '_' participants.participant_id{i}];
    parameters(i).age=participants.age(i);
    parameters(i).sex=participants.sex{i};
    parameters(i).group=participants.group{i};
    parameters(i).site=dset_name;

    % Volumetrics
    lis=[fileList(contains({fileList.name}',t1_string))];
    t1_img=load_nii([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    dimensions=t1_img.original.hdr.dime.pixdim(1:3);
    t1_img=t1_img.img;
    t1_img=reshape(t1_img,[],1);
    t1_sum=sum(t1_img);
    parameters(i).parenchyma=(t1_sum*dimensions(1)*dimensions(2)*dimensions(3));

    lis=[fileList(contains({fileList.name}',gm_nii_string))];
    gm_img=load_nii([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    dimensions=gm_img.original.hdr.dime.pixdim(1:3);
    gm_img=gm_img.img;
    gm_img=reshape(gm_img,[],1);
    gm_sum=sum(gm_img>0.5);
    parameters(i).gm=(gm_sum*dimensions(1)*dimensions(2)*dimensions(3));

    lis=[fileList(contains({fileList.name}',wm_nii_string))];
    wm_img=load_nii([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    dimensions=wm_img.original.hdr.dime.pixdim(1:3);
    wm_img=wm_img.img;
    wm_img=reshape(wm_img,[],1);
    wm_sum=sum(wm_img>0.5);
    parameters(i).wm=(wm_sum*dimensions(1)*dimensions(2)*dimensions(3));

    lis=[fileList(contains({fileList.name}',first_string))];
    first_img=load_nii([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    dimensions=first_img.original.hdr.dime.pixdim(1:3);
    first_img=first_img.img;
    first_img=reshape(first_img,[],1);
    for j=1:length(subcortical_names)
        subcorticals_volume(j,1)=sum(first_img==subcortical_labels(j))*dimensions(1)*dimensions(2)*dimensions(3);
    end
    for j=1:length(subcortical_names)
        parameters(i).(subcortical_names{j})=subcorticals_volume(j);
    end



    % Correlation matrix
    corrmat=atanh(corr(squeeze(ts_clean_reduced(i,:,:))'));
    corrmat(isinf(corrmat))=1;
    parameters(i).corrmat=corrmat;

    % Thickness
    thickness_path_left=strrep(thickness_string,'HEMI','L');
    lis=[fileList(contains({fileList.name}',thickness_path_left))];
    thickness_left=gifti([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    thickness_left=thickness_left.cdata;

    thickness_path_right=strrep(thickness_string,'HEMI','R');
    lis=[fileList(contains({fileList.name}',thickness_path_right))];
    thickness_right=gifti([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    thickness_right=thickness_right.cdata;

    parameters(i).thickness_mp=[thickness_left;thickness_right];

    thickness_roi=zeros(max(glasser_label),1);
    for j=1:max(glasser_label)
        thickness_roi(j,1)=mean(parameters(i).thickness_mp(glasser_label==j),'omitmissing');
    end
    parameters(i).thickness_roi_mp=thickness_roi;

    % Curvature
    curvature_path_left=strrep(curvature_string,'HEMI','L');
    lis=[fileList(contains({fileList.name}',curvature_path_left))];
    curvature_left=gifti([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    curvature_left=curvature_left.cdata;

    curvature_path_right=strrep(curvature_string,'HEMI','R');
    lis=[fileList(contains({fileList.name}',curvature_path_right))];
    curvature_right=gifti([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    curvature_right=curvature_right.cdata;

    parameters(i).curvature_mp=[curvature_left;curvature_right];

    curvature_roi=zeros(max(glasser_label),1);
    for j=1:max(glasser_label)
        curvature_roi(j,1)=mean(parameters(i).curvature_mp(glasser_label==j),'omitmissing');
    end
    parameters(i).curvature_roi_mp=curvature_roi;

    % GD
    lis=[fileList(contains({fileList.name}',gd_string))];
    gd=gifti([lis(contains({lis.folder}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    parameters(i).gd=gd.cdata;

    % Cortical morphology from mri_anatomical_stats
    morphology_string_left=strrep(morphology_string,'HEMI','lh');
    lis=[fileList(contains({fileList.name}',morphology_string_left))];
    morph_left=readmatrix([lis(contains({lis.folder}',participants.participant_id{i})).folder '/' lis(contains({lis.folder}',participants.participant_id{i})).name], 'FileType','text');
    morphology_string_right=strrep(morphology_string,'HEMI','rh');
    lis=[fileList(contains({fileList.name}',morphology_string_right))];
    morph_right=readmatrix([lis(contains({lis.folder}',participants.participant_id{i})).folder '/' lis(contains({lis.folder}',participants.participant_id{i})).name], 'FileType','text');
    morph=[morph_left(2:end,:);morph_right(2:end,:)];
    parameters(i).surface_area=morph(:,2);
    parameters(i).volume=morph(:,3);
    parameters(i).thickness_fs=morph(:,4);
    parameters(i).thickness_fs_sd=morph(:,5);
    parameters(i).rectified_curvature=morph(:,6);
    parameters(i).rectified_curvature_gaussian=morph(:,7);
    parameters(i).folding_index=morph(:,8);
    parameters(i).intrinsic_curvature=morph(:,9);



    if sum(isnan(corrmat(:))) ==0
        % Gradients and their explained variance
        gradfit=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
        gradfit=gradfit.fit(corrmat,'reference',reference_gradients_reduced,'sparsity',80);
        parameters(i).gradients_sp80_na=gradfit.aligned{1};
        parameters(i).lambdas_sp80_na=gradfit.lambda{1}./sum(gradfit.lambda{1});

        gradfit=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
        gradfit=gradfit.fit(corrmat,'reference',reference_gradients_reduced,'sparsity',90);
        parameters(i).gradients_sp90_na=gradfit.aligned{1};
        parameters(i).lambdas_sp90_na=gradfit.lambda{1}./sum(gradfit.lambda{1});

        gradfit=GradientMaps('kernel','cs','approach','dm','alignment','pa','n_components',10);
        gradfit=gradfit.fit(corrmat,'reference',reference_gradients_reduced,'sparsity',90);
        parameters(i).gradients_sp90_cs=gradfit.aligned{1};
        parameters(i).lambdas_sp90_cs=gradfit.lambda{1}./sum(gradfit.lambda{1});
    else
        parameters(i).gradients_sp80_na=zeros(size(reference_gradients_reduced));
        parameters(i).lambdas_sp80_na=zeros(17,1);

        parameters(i).gradients_sp90_na=zeros(size(reference_gradients_reduced));
        parameters(i).lambdas_sp90_na=zeros(17,1);

        parameters(i).gradients_sp90_cs=zeros(size(reference_gradients_reduced));
        parameters(i).lambdas_sp90_cs=zeros(17,1);
    end

end



    %% get the volumetrics old way
    % Volumetrics
    % lis=[fileList(contains({fileList.name}',t1_string))];
    % t1_img=load_nii([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    % dimensions=t1_img.original.hdr.dime.pixdim(1:3);
    % t1_img=t1_img.img;
    % t1_img=reshape(t1_img,[],1);
    % t1_sum=sum(t1_img);
    % parameters(i).parenchyma=(t1_sum*dimensions(1)*dimensions(2)*dimensions(3));
    %
    % lis=[fileList(contains({fileList.name}',gm_nii_string))];
    % gm_img=load_nii([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    % dimensions=gm_img.original.hdr.dime.pixdim(1:3);
    % gm_img=gm_img.img;
    % gm_img=reshape(gm_img,[],1);
    % gm_sum=sum(gm_img>0.5);
    % parameters(i).gm=(gm_sum*dimensions(1)*dimensions(2)*dimensions(3));
    %
    % lis=[fileList(contains({fileList.name}',wm_nii_string))];
    % wm_img=load_nii([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    % dimensions=wm_img.original.hdr.dime.pixdim(1:3);
    % wm_img=wm_img.img;
    % wm_img=reshape(wm_img,[],1);
    % wm_sum=sum(wm_img>0.5);
    % parameters(i).wm=(wm_sum*dimensions(1)*dimensions(2)*dimensions(3));
    %
    % lis=[fileList(contains({fileList.name}',first_string))];
    % first_img=load_nii([lis(contains({lis.name}',participants.participant_id{i})).folder '/' lis(contains({lis.name}',participants.participant_id{i})).name]);
    % dimensions=first_img.original.hdr.dime.pixdim(1:3);
    % first_img=first_img.img;
    % first_img=reshape(first_img,[],1);
    % for j=1:length(subcortical_names)
    %     subcorticals_volume(j,1)=sum(first_img==subcortical_labels(j))*dimensions(1)*dimensions(2)*dimensions(3);
    % end
    % parameters(i).subcorticals=subcorticals_volume;


