
%% I. Convert DICOM --> Nifti
 
% Do everything in parallel
clear all
close all

% access SQL interface
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';

% Load database file.
mksqlite('open', 'ctkDICOM.sql' ) ;
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);
series_all = mksqlite(['select * from Series where SeriesDescription like ''%ADC%'' order by SeriesDate ASC']);
series_DWI_DTI = mksqlite(['select * from Series where '...
        '(SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
        'or SeriesDescription like ''%Fractional Aniso%'' '... % Includes FA from DTI
        'or SeriesDescription like ''%Average DC%'') '...  % Includes ADC from DTI
        'order by SeriesNumber ASC']);

n_studies_all = size(studies_all,1);
n_series_all = size(series_all,1);
n_series_DWI_DTI = size(series_DWI_DTI,1);

fid = fopen([script02_prefix 'DICOM_to_NIfTI.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_series_DWI_DTI 
    fprintf(fid, [' job' num2str(ii)]);
end

zz=1;
for ii = 1:n_studies_all
    
    StudyInstanceUID = studies_all(ii).StudyInstanceUID;
    temp_series = mksqlite(['select * from Series where StudyInstanceUID = ''' StudyInstanceUID ''' '...
            'and (SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
            'or SeriesDescription like ''%Fractional Aniso%'' '... % Includes FA from DTI
            'or SeriesDescription like ''%Average DC%'') '...  % Includes ADC from DTI
            'order by SeriesNumber ASC']);
    
    for jj = 1:numel(temp_series) % will atuomatically pass over if no content in temp_series
        SeriesInstanceUID = temp_series(jj).SeriesInstanceUID;
        SeriesDescription = strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');        
        file_infix = sprintf(['%02d_' SeriesDescription '_' SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID
        fprintf(fid,['\n\njob' num2str(zz) ':\n']);

        images   = mksqlite(['select * from Images where SeriesInstanceUID = ''' SeriesInstanceUID '''' ]);
        pieces = strsplit(images(1).Filename,'/'); % Only works in Matlab 2014a
        halves = strsplit(images(1).Filename,pieces(end)); % Use end fragment to figure out pathname
        pathname = halves{1,1};

        % Convert 1 series to Nifti - Command Source_Directory Output_Filename SeriesInstanceUID(to ensure only that series is used)
        pathname=strrep(pathname,' ','_');
        fprintf(fid, ['\tDicomSeriesReadImageWrite2 ' pathname ' '...
            script02_prefix '00_radpath_raw/radpath_raw_' file_infix '.nii.gz '...
            SeriesInstanceUID '\n']);
        zz=zz+1;
    end
end
    
fclose(fid);
system(['make -j 8 -f ' script02_prefix 'DICOM_to_NIfTI.makefile'])

mksqlite('close' ) ;






%% II. Register ADC map to reference, then apply to DWI-derived maps (ADC, eADC)
%% Register Average DC map to reference, then apply to DTI-Derived maps (Avg DC, FA)
tic
clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';

mksqlite('open', 'ctkDICOM.sql');
series_all = mksqlite(['select * from Series order by SeriesDate ASC']);
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);
series_DWI_DTI = mksqlite(['select * from Series where '...
        '(SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
        'or SeriesDescription like ''%Fractional Aniso%'' '... % Includes FA from DTI
        'or SeriesDescription like ''%Average DC%'') '...  % Includes ADC from DTI
        'order by SeriesNumber ASC']);
    
    
n_series_all = size(series_all,1);
n_studies_all = size(studies_all,1);
n_SWAN = size(series_SWAN,1);
n_T1 = size(series_T1,1);
n_T1post = size(series_T1post,1);
n_T2 = size(series_T2,1);
n_T2star = size(series_T2star,1);
n_T2FLAIR = size(series_T2FLAIR,1);
n_series_DWI_DTI = numel(series_DWI_DTI);

n_register = n_series_DWI_DTI;% Leave out the fixed image

%%%%%%%%%%%%%%%% 1. 
fid = fopen([script02_prefix 'Mass_Register.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_register
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

zz=1;
for ii=1:n_studies_all
    fixed_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%FLAIR%'' ']);
    moving_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%Apparent Diffusion Coefficient%'')'... % Includes ADC and eADC from DWI
        'and SeriesDescription not like ''%Exponential%'' '... % Exclude eADC
        'order by SeriesNumber ASC']);
    
    Fixed_SeriesInstanceUID = fixed_series.SeriesInstanceUID;
    Moving_SeriesInstanceUID = moving_series.SeriesInstanceUID;
    Fixed_SeriesDescription = strrep(strrep(strrep(strrep(strrep(fixed_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    Moving_SeriesDescription = strrep(strrep(strrep(strrep(strrep(moving_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');   
    Fixed_file_infix = sprintf(['%02d_' Fixed_SeriesDescription '_' Fixed_SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID
    Moving_file_infix = sprintf(['%02d_' Moving_SeriesDescription '_' Moving_SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID

    fprintf(fid,['\n\njob' num2str(zz) ':\n']);

    % 1. Perform rough down-sampled registration as beginning step to help
    % initialize later registration
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
        script03_prefix '01_RInput/Affine_Fixed_' Fixed_file_infix '.nii.gz 4 4 4 1\n']); % Down-sample fixed image
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        script03_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz '...
        script03_prefix '01_RInput/Affine_Moving_' Moving_file_infix '.nii.gz 4 4 4 1\n']); %Down-sample moving image
    fprintf(fid,['\t$(ANTSPATH)/antsAffineInitializer 3 '...
        script03_prefix '01_RInput/Affine_Fixed_' Fixed_file_infix '.nii.gz '...
        script03_prefix '01_RInput/Affine_Moving_' Moving_file_infix '.nii.gz '...
        script03_prefix '01_RInput/InitialAffine_' Moving_file_infix '.mat 15 0.1 0 10\n']); % Rough registration

    % 2. Create Laplacian (edge detection) version of images to assist in registration
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script03_prefix '01_RInput/Laplacian_' Fixed_file_infix '.nii.gz '...
        'Laplacian '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz \n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script03_prefix '01_RInput/Laplacian_' Moving_file_infix '.nii.gz '...
        'Laplacian '...
        script03_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz 1.5 1\n']);

    % 3. Perform registration 
    fprintf(fid,['\t$(ANTSPATH)/antsRegistration -d 3 -u 1 -w [0.025, 0.975] '...
        '-o ' script03_prefix '02_ROutput/RegistrationOutput_' Moving_file_infix '_ '...
        '-r ' script03_prefix '01_RInput/InitialAffine_' Moving_file_infix '.mat -z 1 --float 0 '...
        '-m MI[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
              script03_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m MI[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
              script03_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
              '\n']);
%         '-m CC[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
%               script_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz,0.5,4] '...
%         '-m CC[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
%              script_prefix '01_RInput/Laplacian_' Moving_file_infix '.nii.gz,0.5,4] -c [50x10x0,1e-9,15] -t SyN[0.1,3,0] -f 4x2x1 -s 2x1x0'...
%                 '\n']);

    % 4. Apply registration 
    fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
        '-i ' script03_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz '...
        '-o ' script03_prefix '03_Warped/Warped_' Moving_file_infix '.nii.gz '...
        '-r ' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
        '-n NearestNeighbor '... % Interpolation method for resampling
        '-t [' script03_prefix '02_ROutput/RegistrationOutput_' Moving_file_infix '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
        '\n']); 

%             '-t ' script_prefix '02_ROutput/RegistrationOutput_' Moving_file_infix '_1Warp.nii.gz --float 0'... % Deformation Field Transform
%             '\n']); 


    series_DWI_maps = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%exponential%'' '... % eADC map
        'order by SeriesDate, SeriesInstanceUID ASC']);
    n_series_DWI_maps = size(series_DWI_maps,1);

    for jj = 1:n_series_DWI_maps
        Moving_SeriesInstanceUID_maps = series_DWI_maps(jj).SeriesInstanceUID;
        Moving_SeriesDescription_maps = strrep(strrep(strrep(strrep(strrep(series_DWI_maps(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');      
        Moving_file_infix_maps = sprintf(['%02d_' Moving_SeriesDescription_maps '_' Moving_SeriesInstanceUID_maps],ii); % patient #, Descrip, SeriesInstanceUID 
        fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
                '-i ' script03_prefix '00_radpath_raw/' Moving_file_infix_maps '.nii.gz '...
                '-o ' script03_prefix '03_Warped/Warped_' Moving_file_infix_maps '.nii.gz '... 
                '-r ' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
                '-n NearestNeighbor '... % Interpolation method for resampling
                '-t [' script03_prefix '02_ROutput/RegistrationOutput_' Moving_file_infix '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
                '\n']); 
    end
    
    zz=zz+1;
    end
end

for ii=1:n_studies_all
    fixed_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%FLAIR%'' ']);
    moving_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%Average DC%'') '...  % Includes ADC from DTI
        'order by SeriesNumber ASC']);
    
    Fixed_SeriesInstanceUID = fixed_series.SeriesInstanceUID;
    Moving_SeriesInstanceUID = moving_series.SeriesInstanceUID;
    Fixed_SeriesDescription = strrep(strrep(strrep(strrep(strrep(fixed_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    Moving_SeriesDescription = strrep(strrep(strrep(strrep(strrep(moving_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star'); 
    Fixed_file_infix = sprintf(['%02d_' Fixed_SeriesDescription '_' Fixed_SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID
    Moving_file_infix = sprintf(['%02d_' Moving_SeriesDescription '_' Moving_SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID

    fprintf(fid,['\n\njob' num2str(zz) ':\n']);

    % 1. Perform rough down-sampled registration as beginning step to help
    % initialize later registration
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
        script03_prefix '01_RInput/Affine_Fixed_' Fixed_file_infix '.nii.gz 4 4 4 1\n']); % Down-sample fixed image
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        script03_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz '...
        script03_prefix '01_RInput/Affine_Moving_' Moving_file_infix '.nii.gz 4 4 4 1\n']); %Down-sample moving image
    fprintf(fid,['\t$(ANTSPATH)/antsAffineInitializer 3 '...
        script03_prefix '01_RInput/Affine_Fixed_' Fixed_file_infix '.nii.gz '...
        script03_prefix '01_RInput/Affine_Moving_' Moving_file_infix '.nii.gz '...
        script03_prefix '01_RInput/InitialAffine_' Moving_file_infix '.mat 15 0.1 0 10\n']); % Rough registration

    % 2. Create Laplacian (edge detection) version of images to assist in registration
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script03_prefix '01_RInput/Laplacian_' Fixed_file_infix '.nii.gz '...
        'Laplacian '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz \n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script03_prefix '01_RInput/Laplacian_' Moving_file_infix '.nii.gz '...
        'Laplacian '...
        script03_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz 1.5 1\n']);

    % 3. Perform registration 
    fprintf(fid,['\t$(ANTSPATH)/antsRegistration -d 3 -u 1 -w [0.025, 0.975] '...
        '-o ' script03_prefix '02_ROutput/RegistrationOutput_' Moving_file_infix '_ '...
        '-r ' script03_prefix '01_RInput/InitialAffine_' Moving_file_infix '.mat -z 1 --float 0 '...
        '-m MI[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
              script03_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m MI[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
              script03_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
              '\n']);
%         '-m CC[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
%               script_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz,0.5,4] '...
%         '-m CC[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
%              script_prefix '01_RInput/Laplacian_' Moving_file_infix '.nii.gz,0.5,4] -c [50x10x0,1e-9,15] -t SyN[0.1,3,0] -f 4x2x1 -s 2x1x0'...
%                 '\n']);

    % 4. Apply registration 
    fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
        '-i ' script03_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz '...
        '-o ' script03_prefix '03_Warped/Warped_' Moving_file_infix '.nii.gz '...
        '-r ' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
        '-n NearestNeighbor '... % Interpolation method for resampling
        '-t [' script03_prefix '02_ROutput/RegistrationOutput_' Moving_file_infix '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
        '\n']); 

%             '-t ' script_prefix '02_ROutput/RegistrationOutput_' Moving_file_infix '_1Warp.nii.gz --float 0'... % Deformation Field Transform
%             '\n']); 


    series_DTI_maps = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%Fractional Aniso%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);
    n_series_DTI_maps = size(series_DTI_maps,1);

    for jj = 1:n_series_DTI_maps
        Moving_SeriesInstanceUID_maps = series_DTI_maps(jj).SeriesInstanceUID;
        Moving_SeriesDescription_maps = strrep(strrep(strrep(strrep(strrep(strrep(series_DTI_maps(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star'),'.','');
        Moving_file_infix_maps = sprintf(['%02d_' Moving_SeriesDescription_maps '_' Moving_SeriesInstanceUID_maps],ii); % patient #, Descrip, SeriesInstanceUID 
        fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
                '-i ' script03_prefix '00_radpath_raw/' Moving_file_infix_maps '.nii.gz '...
                '-o ' script03_prefix '03_Warped/Warped_' Moving_file_infix_maps '.nii.gz '... 
                '-r ' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
                '-n NearestNeighbor '... % Interpolation method for resampling
                '-t [' script03_prefix '02_ROutput/RegistrationOutput_' Moving_file_infix '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
                '\n']); 
    end
    
    zz=zz+1;

    end
end


fclose(fid);
system(['make -j 8 -f ' script02_prefix 'Mass_Register.makefile'])

time_s = toc
time_min = time_s/60
time_hrs = time_min/60 % 0.6 hours

mksqlite('close' ) ;



%% Check results - make sure image volumes have no hiccups
ii=4 
fixed_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%FLAIR%'' ']);
series_ADC = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and (SeriesDescription like ''%Apparent Diffusion Coefficient%'')'... % Includes ADC and eADC from DWI
    'and SeriesDescription not like ''%Exponential%'' '... % Exclude eADC
    'order by SeriesNumber ASC']);
series_eADC = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and (SeriesDescription like ''%Exponential%'')'... % Includes eADC
    'order by SeriesNumber ASC']);
series_FA = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and (SeriesDescription like ''%Fractional Aniso%'')'... % Includes FA from DTI
    'order by SeriesNumber ASC']); 
series_AvgDC = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and (SeriesDescription like ''%Average DC%'') '...  % Includes ADC from DTI
    'order by SeriesNumber ASC']);     

Fixed_SeriesInstanceUID = fixed_series.SeriesInstanceUID;
  
ADC_SeriesInstanceUID = series_ADC.SeriesInstanceUID; % good, except #3
ADC_SeriesDescription = strrep(strrep(strrep(strrep(strrep(series_ADC.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');  
ADC_file_infix = sprintf(['%02d_' ADC_SeriesDescription '_' ADC_SeriesInstanceUID],ii); %

eADC_SeriesInstanceUID = series_eADC.SeriesInstanceUID; 
eADC_SeriesDescription = strrep(strrep(strrep(strrep(strrep(series_eADC.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
eADC_file_infix = sprintf(['%02d_' eADC_SeriesDescription '_' eADC_SeriesInstanceUID],ii); %

FA_SeriesInstanceUID = series_FA.SeriesInstanceUID; 
FA_SeriesDescription = strrep(strrep(strrep(strrep(strrep(series_FA.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');  
FA_file_infix = sprintf(['%02d_' FA_SeriesDescription '_' FA_SeriesInstanceUID],ii); %

AvgDC_SeriesInstanceUID = series_AvgDC.SeriesInstanceUID; 
AvgDC_SeriesDescription = strrep(strrep(strrep(strrep(strrep(series_AvgDC.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');  
AvgDC_file_infix = sprintf(['%02d_' AvgDC_SeriesDescription '_' AvgDC_SeriesInstanceUID],ii); %

system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script02_prefix '10_Warped/Warped_' ADC_file_infix '.nii.gz '...
    '-o ' script01_prefix '00_radpath_raw/radpath_raw_' Fixed_SeriesInstanceUID '.nii.gz '])




system(['vglrun /opt/apps/SLICER/Slicer-4.3.1-linux-amd64/Slicer'])
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap'])

% View Histogram matched Brain-extracted volumes
system(['vglrun itksnap -g 00_radpath_raw/radpath_raw_' SeriesInstanceUID '.nii.gz'])

