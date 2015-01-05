
 
% Do everything in parallel
clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';
LPBA40_Template_location= '/mnt/data/scratch/igilab/jslin1/LPBA40_Template.nii.gz';
LPBA40_mask_location = '/mnt/data/scratch/igilab/jslin1/LPBA40_mask.nii.gz';

mksqlite('open', 'ctkDICOM.sql' );
series_all = mksqlite(['select * from Series order by SeriesDate ASC']);
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);
series_SWAN = mksqlite(['select * from Series where SeriesDescription like ''%SWAN%'' order by SeriesDate ASC']);%14, none for patient 7
series_T1 = mksqlite(['select * from Series where SeriesDescription like ''%T1%'' and (SeriesDescription not like ''%STEALTH%'' and SeriesDescription not like ''%+C%'') order by SeriesDate ASC']); %17
series_T1post = mksqlite(['select * from Series where SeriesDescription like ''%T1%'' and (SeriesDescription like ''%STEALTH%'' or SeriesDescription like ''%+C%'') order by SeriesDate ASC']); %17
series_T2 = mksqlite(['select * from Series where SeriesDescription like ''%T2%'' and SeriesDescription not like ''%*%'' and SeriesDescription not like ''%FLAIR%'' order by SeriesDate ASC']); %17
series_T2star = mksqlite(['select * from Series where SeriesDescription like ''%T2*%'' order by SeriesDate ASC']); %17
series_T2FLAIR = mksqlite(['select * from Series where SeriesDescription like ''%FLAIR%'' order by SeriesDate ASC']); %17
series_T1_T2_SWAN = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%T1%'' '... % Includes T1, T1post
        'or SeriesDescription like ''%T2%'' '... % Includes T2, T2star, T2FLAIR
        'or SeriesDescription like ''%SWAN%'' '... % Includes SWAN Image volume
        'order by SeriesDate ASC']);
    
prefix = script01_prefix;
series_T2_register = series_T2([2:8 10:size(series_T2,1)]); % Ptno 1, 7 - High rez T2 post
series_T2_HistMatch = series_T2([1 3:7 9:size(series_T2,1)]); % Ptno 1, 7 - lower rez T2

%% I. Convert DICOM --> Nifti
n_cores = 9;
dicom_series = series_T1_T2_SWAN;
makefile = 'DICOM_to_NIfTI.makefile';
fid = fopen([prefix makefile],'w');
fprintf(fid, 'all:');
for ii=1:size(dicom_series,1)
    fprintf(fid, [' job' num2str(ii)]);
end

for jj = 1:size(dicom_series,1) 
    fprintf(fid,['\n\njob' num2str(jj) ':\n']);

    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(dicom_series(jj));
    images   = mksqlite(['select * from Images where SeriesInstanceUID = ''' UID '''' ]);
    pieces = strsplit(images(1).Filename,'/'); % Only works in Matlab 2014a
    halves = strsplit(images(1).Filename,pieces(end)); % Use end fragment to figure out pathname
    pathname = strrep(halves{1,1},' ','_');

    % Convert 1 series to Nifti - Command Source_Directory Output_Filename SeriesInstanceUID(to ensure only that series is used)
    fprintf(fid, ['\tDicomSeriesReadImageWrite2 '...
        pathname ' '...
        prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz '...
        UID '\n']);
end

fclose(fid);
system(['make -j ' num2str(n_cores) ' -f ' prefix makefile])

%% II. Perform Histogram Matching via these steps:
n_cores = 10;


%% 1. N4 correction
tic
N4(series_T1,       prefix, 'T1_N4.makefile',      n_cores)
N4(series_T1post,   prefix, 'T1post_N4.makefile',  n_cores)
N4(series_T2,       prefix, 'T2_N4.makefile',      n_cores)
N4(series_T2star,   prefix, 'T2star_N4.makefile',  n_cores)
N4(series_T2FLAIR,  prefix, 'T2FLAIR_N4.makefile', n_cores)
N4(series_SWAN,     prefix, 'SWAN_N4.makefile',    n_cores)


%% 2. Registration to T2 Volume
% Re-do for ptno's [1 7] 
tic
fixed_series = series_T2_register; % Exclude Ax T2 24 slices sequence for ptno 7
Register1(prefix, 'T1_Register.makefile',      n_cores, fixed_series, series_T1,      series_T1)
Register1(prefix, 'T1post_Register.makefile',  n_cores, fixed_series, series_T1post,  series_T1post)
Register1(prefix, 'T2star_Register.makefile',  n_cores, fixed_series, series_T2star,  series_T2star)
Register1(prefix, 'T2FLAIR_Register.makefile', n_cores, fixed_series, series_T2FLAIR, series_T2FLAIR)
Register1(prefix, 'SWAN_Register.makefile',    n_cores, fixed_series, series_SWAN,    series_SWAN)


%% 3. Create Masks
MD_parameter = 5;
MC_parameter = 4;
Create_Mask(prefix, 'T2_Create_Mask.makefile'      , n_cores, series_T2,      MD_parameter, MC_parameter)
time_s = toc
time_min = time_s/60
time_hrs = time_min/60 % 27.6 hrs



% Do I need these?
Create_Mask(prefix, 'T1_Create_Mask.makefile'      , n_cores, series_T1,      MD_parameter, MC_parameter)
Create_Mask(prefix, 'T1post_Create_Mask.makefile'  , n_cores, series_T1post,  MD_parameter, MC_parameter)
Create_Mask(prefix, 'T2star_Create_Mask.makefile'  , n_cores, series_T2star,  MD_parameter, MC_parameter)
tic
Create_Mask(prefix, 'T2FLAIR_Create_Mask.makefile' , n_cores, series_T2FLAIR, MD_parameter, MC_parameter)
time_hrs = toc/60/60 
Create_Mask(prefix, 'SWAN_Create_Mask.makefile'    , n_cores, series_SWAN,    MD_parameter, MC_parameter)


% Check mask - how well it matches head volume    
ii=4;
[ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_T1(ii));
system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz '...
    '-s ' script01_prefix '04_Masks/Atropos_' ptno_Descrip_UID '.nii.gz '])


% 4. Extract brain, Truncate Intensity
Extract_Truncate(prefix, fixed_series, series_T1)
Extract_Truncate(prefix, fixed_series, series_T1post)
Extract_Truncate(prefix, fixed_series, series_T2)
Extract_Truncate(prefix, fixed_series, series_T2star)
Extract_Truncate(prefix, fixed_series, series_T2FLAIR)
Extract_Truncate(prefix, fixed_series, series_SWAN)



%% 4. Generate Cumulative Histogram
[T1_cume] = Generate_Cume(prefix, series_T1);
save_nii(T1_cume,[prefix 'T1_cume.nii.gz'])
clear T1_cume

[T1post_cume] = Generate_Cume(prefix, series_T1post);
save_nii(T1post_cume,[prefix 'T1post_cume.nii.gz'])
clear T1post_cume 

[T2_cume] = Generate_Cume(prefix, series_T2_HistMatch); % Only use the ptno 7 AxT2 24 slice
save_nii(T2_cume,[prefix 'T2_cume.nii.gz')
clear T2_cume 

[T2star_cume] = Generate_Cume(prefix, series_T2star);
save_nii(T2star_cume,[prefix  'T2star_cume.nii.gz'])
clear T2star_cume 

[T2FLAIR_cume] = Generate_Cume(prefix, series_T2FLAIR);
save_nii(T2FLAIR_cume,[prefix 'T2FLAIR_cume.nii.gz'])
clear T2FLAIR_cume 

[SWAN_cume] = Generate_Cume(prefix, series_SWAN);
save_nii(SWAN_cume, [script01_prefix 'SWAN_cume.nii.gz'])
clear SWAN_cume



%% 5. Perform Histogram Matching
n_intensity_levels = 255;
n_landmarks = 64;
HistMatch(prefix, 'T1_HistMatch.makefile',      'T1_cume.nii.gz',      n_cores, series_T1,      n_intensity_levels, n_landmarks)
HistMatch(prefix, 'T1post_HistMatch.makefile',  'T1post_cume.nii.gz',  n_cores, series_T1post,  n_intensity_levels, n_landmarks)
HistMatch(prefix, 'T2_HistMatch.makefile',      'T2_cume.nii.gz',      n_cores, series_T2,      n_intensity_levels, n_landmarks)
HistMatch(prefix, 'T2star_HistMatch.makefile',  'T2star_cume.nii.gz',  n_cores, series_T2star,  n_intensity_levels, n_landmarks)
HistMatch(prefix, 'T2FLAIR_HistMatch.makefile', 'T2FLAIR_cume.nii.gz', n_cores, series_T2FLAIR, n_intensity_levels, n_landmarks)
HistMatch(prefix, 'SWAN_HistMatch.makefile',    'SWAN_cume.nii.gz',    n_cores, series_SWAN,    n_intensity_levels, n_landmarks)


% Check results - make sure image volumes have no hiccups

% Generate histogram 1
img = load_nii([script01_prefix '00_radpath_raw/radpath_raw_' UID '.nii.gz']); 
img_values = double(img.img(img.img(:)>0));
img_mean = mean(img_values,1);
img_std = std(img_values,0,1);
h = figure;
hist(img_values)
title([ Descrip ', mean = ' num2str(img_mean) ' +/- ' num2str(img_std)])
saveas(h, [script01 '_Histogram_' Descrip ],'png')

% Generate histogram 2
img = load_nii([script01_prefix '07_HistMatch/Brain_HistMatch_' UID '.nii.gz']); 
img_values = double(img.img(img.img(:)>0));
img_mean = mean(img_values,1);
img_std = std(img_values,0,1);
h = figure;
hist(img_values)
title([ Descrip ', mean = ' num2str(img_mean) ' +/- ' num2str(img_std)])
saveas(h, [script01 '_Histogram_' Descrip ],'png')

% Cumulative Histogram
% Descrip_cume = 'T1_cume';
% Descrip_cume = 'T1post_cume';
% Descrip_cume = 'T2_cume';
% Descrip_cume = 'T2star_cume';
% Descrip_cume = 'T2FLAIR_cume';
% Descrip_cume = 'SWAN_cume';
img = load_nii([script01_prefix Descrip_cume '.nii.gz']); 
img_values = double(img.img(img.img(:)>0));
img_mean = mean(img_values,1);
img_std = std(img_values,0,1);
h = figure;
hist(img_values)
title([ temp_series(jj).SeriesDescription ', mean = ' num2str(img_mean) ' +/- ' num2str(img_std)])
saveas(h, [script01 'Histogram_' Descrip_cume ],'png')

    