
% vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap
% vglrun /opt/apps/amira560/bin/start
%% I. Load stuff
clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql

%% Select Image Volume
ii=1;
prefix = script01_prefix;
moving_series = series_T1(ii);
moving_series = series_T1post(ii);
moving_series = series_T2star(ii);
moving_series = series_T2FLAIR(ii);
moving_series = series_SWAN(ii);
[m_ptno m_Descrip m_UID m_ptno_Descrip_UID] = Generate_Label(moving_series);
conv_factor = 1;
fixed_series = series_T2_register;

ii=1;
prefix = script02_prefix;
moving_series = series_ADC(ii);;
% moving_series = series_FA(ii);
% moving_series = series_AvgDC(ii);
[m_ptno m_Descrip m_UID m_ptno_Descrip_UID] = Generate_Label(moving_series);
conv_factor = 1e-6;
fixed_series = series_T2_register;

ii=1;
prefix = script03_prefix;
moving_series = series_DCE(ii);
moving_series = series_DSC(ii);
moving_series = series_rBV(ii);
[m_ptno m_Descrip m_UID m_ptno_Descrip_UID] = Generate_Label(moving_series);
fixed_series = series_T2FLAIR;
images = mksqlite(['select * from Images where SeriesInstanceUID = ''' m_UID ''' order by Filename ASC' ]);
tmp2=dicominfo(images(1).Filename,'dictionary','gems-dicom-dict.txt');
conv_factor = tmp2.Private_0077_1001;

ii=1;
ii=(ii*4)-3;
% moving_series = series_Ktrans(ii);
% moving_series = series_Vp(ii);
% moving_series = series_AUC(ii);
[m_ptno m_Descrip m_UID m_ptno_Descrip_UID] = Generate_Label(moving_series);
fixed_series = series_T2FLAIR;
images = mksqlite(['select * from Images where SeriesInstanceUID = ''' m_UID ''' order by Filename ASC' ]);
tmp2=dicominfo(images(1).Filename,'dictionary','gems-dicom-dict.txt');
conv_factor = tmp2.Private_0077_1001;

%% Mask check 
system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' prefix '01_N4/N4_' m_ptno_Descrip_UID '.nii.gz '...
    '-s ' prefix '04_Masks/Mask_' m_ptno_Descrip_UID '.nii.gz '])

system(['vglrun /opt/apps/amira560/bin/start'])
%% Registration Comparison

[f_field f_study]= find(cell2mat(cellfun(@(x) isequal(x,moving_series.StudyInstanceUID),struct2cell(fixed_series),'UniformOutput',false)));
[f_ptno f_Descrip f_UID f_ptno_Descrip_UID] = Generate_Label(fixed_series(f_study));

% Registration comparison
f_title = [script01_prefix '00_radpath_raw/radpath_raw_' f_ptno_Descrip_UID '.nii.gz'];
f_title = [script01_prefix '01_N4/N4_' f_ptno_Descrip_UID '.nii.gz'];
f_title = [script01_prefix '04_Warped/Warped_' f_ptno_Descrip_UID '.nii.gz'];
m_title = [prefix '04_Warped/Warped_' m_ptno_Descrip_UID '.nii.gz'];
system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' f_title ' '...
    '-o ' m_title])

system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap'])



%% Histogram 
m_title = ['00_radpath_raw/radpath_raw_' m_ptno_Descrip_UID '.nii.gz'];
m_title = ['06_GaussNorm/GN_' m_ptno_Descrip_UID '.nii.gz'];
save_dir = [prefix 'Histograms/'];
if ~isdir(save_dir)
    mkdir(save_dir)
end

% Matlab
img = load_nii([prefix m_title]); 
% img_values = double(img.img(img.img(:)>0))*conv_factor;
% img_values = double(img.img(img.img(:)>0));
img_values = double(img.img(:));

img_mean = mean(img_values(:),1);
img_std = std(img_values(:),0,1);
h = figure;
% hist(img_values(:),[0:5:5000])
% hist(img_values(:))
xlim([-100 3000])

hist(img_values(:),[0:.1:40])
xlim([-1 25])
title({[ m_ptno ' ' m_Descrip ];[ ' mean = ' num2str(img_mean) ' +/- ' num2str(img_std)]},'Interpreter','none')
xlabel('Intensity Values')
ylabel('# Voxels')

saveas(h, [save_dir m_ptno '_Histogram_' m_Descrip ],'png')

close all

% ITKSnap
system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' prefix m_title])

