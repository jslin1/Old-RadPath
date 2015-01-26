
clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql
lot=3;
prefix = script03_prefix;
n_cores = 8;

%% 000. Copy/Re-number DICOMs for Nifti Conversion
DCE_DSC_Copy_Renumber_DICOMs

%% 00. DICOM --> Nifti
n_cores = 9;
makefile = 'DICOM_to_NIfTI.makefile';
series_sql = series_DCE_DSC_plusmaps;
DICOM_to_Nifti(lot,prefix, makefile, series_sql, n_cores)

%% 0. Fix headers of DCE/DSC maps - spatial location info = applying registration works
DCE_DSC_Fix_Headers

%% 1. N4 correction
N4(series_DCE,       prefix, 'DCE_N4.makefile',      n_cores)
N4(series_DSC,       prefix, 'DSC_N4.makefile',      n_cores)

%% 2. Create Masks
MD_parameter = 5;
MC_parameter = 4;
Create_Mask(lot, prefix, 'DCE_Create_Mask.makefile', n_cores, series_DCE, MD_parameter, MC_parameter)
Create_Mask(lot, prefix, 'DSC_Create_Mask.makefile', n_cores, series_DSC, MD_parameter, MC_parameter)

% vglrun /opt/apps/amira560/bin/start

%% 3. Extract
tic
BO(lot, prefix, series_DCE, series_DCE_plusmaps) 
BO(lot, prefix, series_DSC, series_DSC_plusmaps) 
toc/60 % 8.77 min

%% 4. Register DCE/DSC lasttimept to Reference, then apply to all derived maps
tic
fixed_series = series_T2_register;
Register(lot,prefix, 'DCE_register.makefile', n_cores, fixed_series, series_DCE, series_DCE, series_DCE_plusmaps)
Register(lot,prefix, 'DSC_register.makefile', n_cores, fixed_series, series_DSC, series_DSC, series_DSC_plusmaps)
t_min = toc/60
t_hrs = t_min/60 % 3.94 hrs

% Skull, DCE/DSC = success? it mainly matches the skull?
% Brain --> Brain, DCE/DSC = success
% Brain --> Brain, AUC/rBV = okay, worse than DCE/DSC

% Brain --> Skull, AUC and rBV fail
% Brain --> Skull, DCE and DSC fail

%% 5. Gaussian Normalization
tic
GaussNorm(lot, prefix, series_DCE_plusmaps)
GaussNorm(lot, prefix, series_DSC_plusmaps)
toc/60 % 8.78 min