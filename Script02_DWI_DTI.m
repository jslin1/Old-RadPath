


clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';
mksqlite('open', 'ctkDICOM.sql' );

series_all = mksqlite(['select * from Series order by SeriesDate ASC']);
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);

series_maps = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
        'or SeriesDescription like ''%Fractional Aniso%'' '... % Includes FA from DTI
        'or SeriesDescription like ''%Average DC%'' '...  % Includes ADC from DTI
        'order by SeriesDate ASC']);

series_T2 = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%T2%'' '... % Includes T2, T2star, T2FLAIR
        'and SeriesDescription not like ''%*%'' '... % Exclude T2star
        'and SeriesDescription not like ''%FLAIR%'' '... % Exclude T2 FLAIR
        'order by SeriesDate ASC']);
series_ADC = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
        'and SeriesDescription not like ''%Exponential%'' '... % Exclude eADC
        'order by SeriesDate ASC']);
series_ADC_eADC = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
        'order by SeriesDate ASC']);
series_AvgDC = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%Average DC%'' '...  % Includes ADC from DTI
        'order by SeriesDate ASC']);
series_AvgDC_FA = mksqlite(['select * from Series where '...
       'SeriesDescription like ''%Average DC%'' '...  % Includes ADC from DTI
       'or SeriesDescription like ''%Fractional Aniso%'' '...
       'order by SeriesDate, SeriesInstanceUID ASC']);

prefix = script02_prefix;
series_T2_register = series_T2([2:8 10:size(series_T2,1)]); % T2+C

%% I. DICOM --> Nifti
n_cores = 9;
makefile = 'DICOM_to_NIfTI.makefile';
fid = fopen([prefix makefile],'w');
fprintf(fid, 'all:');
for ii=1:size(series_maps,1)
    fprintf(fid, [' job' num2str(ii)]);
end

for jj = 1:size(series_maps,1) 
    fprintf(fid,['\n\njob' num2str(jj) ':\n']);
    
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_maps(jj));
    images   = mksqlite(['select * from Images where SeriesInstanceUID = ''' UID '''' ]);
    pieces = strsplit(images(1).Filename,'/'); % Only works in Matlab 2014a
    halves = strsplit(images(1).Filename,pieces(end)); % Use end fragment to figure out pathname
    pathname = strrep(halves{1,1},' ','_');

    % Convert 1 series to Nifti - Command Source_Directory Output_Filename SeriesInstanceUID(to ensure only that series is used)
    fprintf(fid, ['\tDicomSeriesReadImageWrite2 ' ...
        pathname ' '...
        prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz '...
        UID '\n']);
end

fclose(fid);
system(['make -j ' num2str(n_cores) ' -f ' prefix makefile])

%% II. Register ADC map to reference, then apply to DWI-derived maps (ADC, eADC)
tic
n_cores = 10;
fixed_series = series_T2_register;
Register2(prefix, 'DWI_register.makefile', n_cores, fixed_series, series_ADC, series_ADC_eADC)
Register2(prefix, 'DTI_register.makefile', n_cores, fixed_series, series_AvgDC, series_AvgDC_FA)
time_s = toc
time_min = time_s/60
time_hrs = time_min/60 % 2.91 hours




%% Check results - make sure image volumes have no hiccups
ii=4
fixed_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%T2%'' '... % Includes T2, T2star, T2FLAIR
        'and SeriesDescription not like ''%*%'' '... % Exclude T2star
        'and SeriesDescription not like ''%FLAIR%'' '... % Exclude T2 FLAIR
        'order by SeriesNumber ASC']);
    
moving_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and (SeriesDescription like ''%Apparent Diffusion Coefficient%'')'... % Includes ADC and eADC from DWI
    'and SeriesDescription not like ''%Exponential%'' '... % Exclude eADC
    'order by SeriesNumber ASC']);
moving_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and (SeriesDescription like ''%Exponential%'')'... % Includes eADC
    'order by SeriesNumber ASC']);
moving_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and (SeriesDescription like ''%Fractional Aniso%'')'... % Includes FA from DTI
    'order by SeriesNumber ASC']); 
moving_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and (SeriesDescription like ''%Average DC%'') '...  % Includes ADC from DTI
    'order by SeriesNumber ASC']);     

[Fixed_ptno Fixed_Descrip Fixed_UID Fixed_ptno_Descrip_UID] = Generate_Label(fixed_series);
[Moving_ptno Moving_Descrip Moving_UID Moving_ptno_Descrip_UID] = Generate_Label(moving_series);

system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '01_N4/Brain_N4_' Fixed_ptno_Descrip_UID '.nii.gz '...
    '-o ' script02_prefix '10_Warped/Warped_' Moving_ptno_Descrip_UID '.nii.gz '])


