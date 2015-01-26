
clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql
lot=2;
prefix = script02_prefix;

%% 0. DICOM --> Nifti
n_cores = 9;
makefile = 'DICOM_to_NIfTI.makefile';
series_sql = series_diffmaps_alone;
DICOM_to_Nifti(lot, prefix, makefile, series_sql, n_cores)

%% 1. Register ADC map to reference, then apply to DWI-derived maps (ADC, eADC)
tic
lot=2;
n_cores = 10;
fixed_series = series_T2_register;
Register(lot, prefix, 'DWI_register.makefile', n_cores, fixed_series, series_ADC, seriesADC, series_ADC_eADC)
Register(lot, prefix, 'DTI_register.makefile', n_cores, fixed_series, series_FA, series_FA, series_AvgDC_FA)
time_min = toc/60
time_hrs = time_min/60 % 2.91 hours

%% 2. Extract Brain
BO(lot, prefix, fixed_series, series_ADC)
BO(lot, prefix, fixed_series, series_eADC)
BO(lot, prefix, fixed_series, series_AvgDC)
BO(lot, prefix, fixed_series, series_FA)

%% 3. Perform Gaussian Normalization
GaussNorm(lot, prefix, series_ADC)
GaussNorm(lot, prefix, series_eADC)
GaussNorm(lot, prefix, series_AvgDC)
GaussNorm(lot, prefix, series_FA)