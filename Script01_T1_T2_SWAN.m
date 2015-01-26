
clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql
lot=1;
prefix = script01_prefix;
n_cores = 8;
fixed_series = series_T2_register; % Exclude Ax T2 24 slices sequence for ptno 7

%% 0. Convert DICOM --> Nifti
n_cores = 9;
series_sql = series_T1post(end)%series_T1_T2_SWAN;
makefile = 'DICOM_to_NIfTI.makefile';
DICOM_to_Nifti(lot, prefix, makefile, series_sql, n_cores)

%% Perform Histogram Matching via these steps:
%% 1. N4 correction
n_cores = 8;
N4(series_T1,       prefix, 'T1_N4.makefile',      n_cores)
N4(series_T1post,   prefix, 'T1post_N4.makefile',  n_cores)
N4(series_T2,       prefix, 'T2_N4.makefile',      n_cores)
N4(series_T2star,   prefix, 'T2star_N4.makefile',  n_cores)
N4(series_T2FLAIR,  prefix, 'T2FLAIR_N4.makefile', n_cores)
N4(series_SWAN,     prefix, 'SWAN_N4.makefile',    n_cores)

%% 2. Create T2, T2FLAIR Masks
n_cores = 8;
MD_parameter = 5;
MC_parameter = 4;
Create_Mask(0, prefix, 'T2_Create_Mask.makefile'      , n_cores, series_T2,      MD_parameter, MC_parameter)
Create_Mask(1, prefix, 'T2FLAIR_Create_Mask.makefile' , n_cores, series_T2FLAIR, MD_parameter, MC_parameter)
t_min = toc/60
t_hrs = toc/60/60 % 27.6 hrs

% vglrun /opt/apps/amira560/bin/start - +C = -C must fix manually

% Do I need these?
% Create_Mask(type,prefix, 'T1_Create_Mask.makefile'      , n_cores, series_T1,      MD_parameter, MC_parameter)
% Create_Mask(type,prefix, 'T1post_Create_Mask.makefile'  , n_cores, series_T1post,  MD_parameter, MC_parameter)
% Create_Mask(type,prefix, 'T2star_Create_Mask.makefile'  , n_cores, series_T2star,  MD_parameter, MC_parameter)
% Create_Mask(type,prefix, 'SWAN_Create_Mask.makefile'    , n_cores, series_SWAN,    MD_parameter, MC_parameter)

%% 3. Registration to T2 Volume
n_cores = 10;
tic
Register(lot, prefix, 'T2_Register.makefile',  n_cores, fixed_series, series_T2_HistMatch([1 7]),  series_T2_HistMatch([1 7]), series_T2_HistMatch([1 7]))
Register(lot, prefix, 'T1_Register.makefile',      n_cores, fixed_series, series_T1,      series_T1,        series_T1)


Register(lot, prefix, 'T1post_Register.makefile',  n_cores, fixed_series, series_T1post,  series_T1post,    series_T1post)
Register(lot, prefix, 'T2star_Register.makefile',  n_cores, fixed_series, series_T2star,  series_T2star,    series_T2star)
Register(lot, prefix, 'T2FLAIR_Register.makefile', n_cores, fixed_series, series_T2FLAIR, series_T2FLAIR,   series_T2FLAIR)
Register(lot, prefix, 'SWAN_Register.makefile',    n_cores, fixed_series, series_SWAN,    series_SWAN,      series_SWAN)
toc/60
toc/60/60

%% 4. Extract brain, Truncate Intensity
tic
BO(lot, prefix, fixed_series, series_T1)
BO(lot, prefix, fixed_series, series_T1post)
BO(0,   prefix, fixed_series, series_T2_register)
BO(lot, prefix, fixed_series, series_T2star)
BO(lot, prefix, fixed_series, series_T2FLAIR)
BO(lot, prefix, fixed_series, series_SWAN)
toc/60

%% 5. Perform Gaussian Normalization
tic
GaussNorm(lot, prefix, series_T1)
GaussNorm(lot, prefix, series_T1post)
GaussNorm(0, prefix, series_T2)
GaussNorm(lot, prefix, series_T2star)
GaussNorm(lot, prefix, series_T2FLAIR)
GaussNorm(lot, prefix, series_SWAN)
toc/60

    