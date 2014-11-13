%% 0a. Extract last timept
clear all
close all

% access SQL interface
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';

% Load database file.
mksqlite('open', 'ctkDICOM.sql' ) ;
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);
series_all = mksqlite(['select * from Series order by SeriesDate, SeriesInstanceUID ASC']);
series_DCE_DSC = mksqlite(['select * from Series '...
    'where SeriesDescription like ''%DCE%'' '...
    'or SeriesDescription like ''%DSC%'' '...
    'order by SeriesDate, SeriesInstanceUID ASC']);

n_studies_all = numel(studies_all);
n_series_all = numel(series_all);
n_series_DCE_DSC = numel(series_DCE_DSC);


for kk = 4%1:n_studies_all
    StudyInstanceUID = studies_all(kk).StudyInstanceUID;
    temp_series = mksqlite(['select * from Series where StudyInstanceUID = ''' StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%DCE%'' '...
        'or SeriesDescription like ''%DSC%'')'...
        'order by SeriesDate, SeriesInstanceUID ASC']);
    n_temp_series = numel(temp_series);

    for jj=1:n_temp_series % for NICE/NNL, Olea DROs
        SeriesInstanceUID = temp_series(jj).SeriesInstanceUID;
        SeriesDescription = strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    
        images   = mksqlite(['select * from Images where SeriesInstanceUID = ''' SeriesInstanceUID ''' order by SOPInstanceUID ASC']);
        pieces = strsplit(images(1).Filename,'/'); % Only works in Matlab 2014a
        halves = strsplit(images(1).Filename,pieces(end)); % Use end fragment to figure out pathname
        pathname1 = halves{1,1};
        pathname2 = [script03_prefix '00_radpath_raw/DICOMs_Re-numbered/' sprintf(['%02d'],kk) '/' sprintf(['%02d'],kk) '_' SeriesDescription '_lasttimept_DICOMs' ];

        if size(dir([pathname1 '*.dcm']),1)>0
            tmp_pathname = dir([pathname1 '*.dcm']);
        elseif size(dir([pathname1 '*IM*']),1)>0
            tmp_pathname = dir([pathname1 '*IM*']);
        end

        file_names = strvcat(tmp_pathname(1:size(tmp_pathname,1)).name); % Generates file names for every file in that directory
        nimages = size(file_names,1); % # of files in that directory

        images(1).Filename
        tmp2=dicominfo(images(1).Filename,'dictionary','gems-dicom-dict.txt');
        LocationsInAcquisition=tmp2.LocationsInAcquisition; % # Slices

        clear TemporalPositionIdentifier  NumberOfTemporalPositions
        disp(['Loading Images'])
        for ii = nimages-200:nimages; % Reads in dicom files - the actual data, not just the filenames
            tmp2=dicominfo(images(ii).Filename,'dictionary','gems-dicom-dict.txt');
            TemporalPositionIdentifier(ii,1)=tmp2.TemporalPositionIdentifier;%?
            NumberOfTemporalPositions(ii,1)=tmp2.NumberOfTemporalPositions; % #timepoints
            if NumberOfTemporalPositions(ii,1)==TemporalPositionIdentifier(ii,1) % Only for last timept
                [ii tmp2.SliceLocation tmp2.InstanceNumber]
                SliceLocation_lasttimept(ii,1)=tmp2.SliceLocation;
                InstanceNumber_lasttimept(ii,1)=tmp2.InstanceNumber;
                if ~isdir([pathname2])
                  mkdir([pathname2]) % Create root image directory
                end
                status = copyfile(images(ii).Filename, [pathname2 '/' sprintf('%05d',tmp2.InstanceNumber)],'f'); % finish this command
            end
        end
        disp(['Done!'])
    end
end

%% 0b. Re-number maps
clear all
close all

% access SQL interface
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';

% Load database file.
mksqlite('open', 'ctkDICOM.sql' ) ;
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);
series_all = mksqlite(['select * from Series order by SeriesDate, SeriesInstanceUID ASC']);
n_studies_all = numel(studies_all);
n_series_all = numel(series_all);


for kk = 4%1:n_studies_all
    temp_series = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(kk).StudyInstanceUID ''' '...
            'and (SeriesDescription like ''%rBV%'' '...
            'or SeriesDescription like ''%rBF%'' '...
            'or SeriesDescription like ''%MTT%'' '...
            'or SeriesDescription like ''%Delay%'' '...
            'or SeriesDescription like ''%Leakage (K2)%'' '...
            'or SeriesDescription like ''%K12%'' '...
            'or SeriesDescription like ''%K21%'' '...
            'or SeriesDescription like ''%Vp%'' '...
            'or SeriesDescription like ''%Ve (distribution volume)%'' '...
            'or SeriesDescription like ''%Wash-in%'' '...
            'or SeriesDescription like ''%Wash-out%'' '...
            'or SeriesDescription like ''%Time to peak%'' '...
            'or SeriesDescription like ''%Area under curve%'' '...
            'or SeriesDescription like ''%Peak enhancement%'')'...
            'order by SeriesDate, SeriesInstanceUID ASC']);
    n_temp_series = numel(temp_series);

    for jj=1:n_temp_series
        SeriesInstanceUID = temp_series(jj).SeriesInstanceUID;
        SeriesDescription = strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
        file_infix = sprintf(['%02d_' SeriesDescription '_' SeriesInstanceUID],kk);

        images   = mksqlite(['select * from Images where SeriesInstanceUID = ''' SeriesInstanceUID ''' order by Filename ASC' ]);
        pieces = strsplit(images(1).Filename,'/'); % Only works in Matlab 2014a
        halves = strsplit(images(1).Filename,pieces(end)); % Use end fragment to figure out pathname
        pathname1 = halves{1,1};
        pathname2 = [script03_prefix '00_radpath_raw/DICOMs_Re-numbered/' sprintf(['%02d'],kk) '/' sprintf(['%02d'],kk) '_' SeriesDescription '_DICOMs' ];

        if size(dir([pathname1 '*.dcm']),1)>0
            tmp_pathname = dir([pathname1 '*.dcm']);
        elseif size(dir([pathname1 '*IM*']),1)>0
            tmp_pathname = dir([pathname1 '*IM*']);
        end

        file_names = strvcat(tmp_pathname(1:size(tmp_pathname,1)).name); % Generates file names for every file in that directory
        nimages = size(file_names,1); % # of files in that directory

        for ii = 1:nimages; % Reads in dicom files - the actual data, not just the filenames
            tmp2=dicominfo(images(ii).Filename,'dictionary','gems-dicom-dict.txt');
            [ii tmp2.SliceLocation tmp2.InstanceNumber]
            SliceLocation_lasttimept(ii,1)=tmp2.SliceLocation;
            InstanceNumber_lasttimept(ii,1)=tmp2.InstanceNumber;
            if ~isdir([pathname2])
              mkdir([pathname2]) % Create root image directory
            end
            status = copyfile(images(ii).Filename, [pathname2 '/' sprintf('%05d',tmp2.InstanceNumber)],'f'); % finish this command
        end
    end
end




%% I. DICOM --> Nifti
% Start by generating database by setting Local Database directory on Slicer's DICOM
% database tool, then importing DICOMs

clear all
close all

% access SQL interface
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';

% Load database file.
mksqlite('open', 'ctkDICOM.sql');
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);
series_all = mksqlite(['select * from Series where SeriesDescription like ''%%'' or SeriesDescription like ''%%'' order by SeriesDate ASC']);
n_studies_all = size(studies_all,1);
n_series_all = size(series_all,1);

fid = fopen([script03_prefix 'DICOM_to_NIfTI.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_series_all 
    fprintf(fid, [' job' num2str(ii)]);
end

zz=1;
for kk = 4;%1:n_studies
    StudyInstanceUID = studies_all(kk).StudyInstanceUID;
    temp_series = mksqlite(['select * from Series where StudyInstanceUID = ''' StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%DCE%'' '...
        'or SeriesDescription like ''%DSC%'' '...
        'or SeriesDescription like ''%rBV%'' '...
        'or SeriesDescription like ''%rBF%'' '...
        'or SeriesDescription like ''%MTT%'' '...
        'or SeriesDescription like ''%Delay%'' '...
        'or SeriesDescription like ''%Leakage (K2)%'' '...
        'or SeriesDescription like ''%K12%'' '...
        'or SeriesDescription like ''%K21%'' '...
        'or SeriesDescription like ''%Vp%'' '...
        'or SeriesDescription like ''%Ve (distribution volume)%'' '...
        'or SeriesDescription like ''%Wash-in%'' '...
        'or SeriesDescription like ''%Wash-out%'' '...
        'or SeriesDescription like ''%Time to peak%'' '...
        'or SeriesDescription like ''%Area under curve%'' '...
        'or SeriesDescription like ''%Peak enhancement%'')'...
        'order by SeriesDate, SeriesInstanceUID ASC']);
    n_temp_series = numel(temp_series);
    
    for jj = 1:n_temp_series
        SeriesInstanceUID = temp_series(jj).SeriesInstanceUID;
        SeriesDescription = strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
        fprintf(fid,['\n\njob' num2str(zz) ':\n']);

        if numel(strfind(SeriesDescription,'DCE'))>0 || numel(strfind(SeriesDescription,'DSC'))>0 % If DCE or DSC found, execute this code
            pathname = [script03_prefix '00_radpath_raw/DICOMs_Re-numbered/' sprintf(['%02d'],kk) '/' sprintf(['%02d'],kk) '_' SeriesDescription '_lasttimept_DICOMs' ];
            file_infix = sprintf(['%02d_' SeriesDescription '_lasttimept_' SeriesInstanceUID],kk);
            fprintf(fid, ['\tDicomSeriesReadImageWrite2 '...
            pathname ' '...
            script03_prefix '00_radpath_raw/radpath_raw_' file_infix '.nii.gz '...
            SeriesInstanceUID '\n']);
            zz=zz+1;
        else
            pathname = [script03_prefix '00_radpath_raw/DICOMs_Re-numbered/' sprintf(['%02d'],kk) '/' sprintf(['%02d'],kk) '_' SeriesDescription '_DICOMs' ];    
            file_infix = sprintf(['%02d_' SeriesDescription '_' SeriesInstanceUID],kk);
            fprintf(fid, ['\tDicomSeriesReadImageWrite2 '...
            pathname ' '...
            script03_prefix '00_radpath_raw/radpath_raw_' file_infix '.nii.gz '...
            SeriesInstanceUID '\n']);
            zz=zz+1;
        end
    end

end

fclose(fid);
system(['make -j 9 -f ' script03_prefix 'DICOM_to_NIfTI.makefile'])

mksqlite('close' ) ;


%% Ia. Fix headers of DCE maps and DSC maps - not that time consuming
ii=4;
series_DCE = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%DCE%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);
series_DCE_maps =     mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%K12%'' '...
        'or SeriesDescription like ''%K21%'' '...
        'or SeriesDescription like ''%Vp%'' '...
        'or SeriesDescription like ''%Ve (distribution volume)%'' '...
        'or SeriesDescription like ''%Wash-in%'' '...
        'or SeriesDescription like ''%Wash-out%'' '...
        'or SeriesDescription like ''%Time to peak%'' '...
        'or SeriesDescription like ''%Area under curve%'' '...
        'or SeriesDescription like ''%Peak enhancement%'')'...
        'order by SeriesDate, SeriesInstanceUID ASC']);
SeriesInstanceUID = series_DCE.SeriesInstanceUID;
SeriesDescription = strrep(strrep(strrep(strrep(strrep(series_DCE.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
DCE_file_infix = sprintf(['%02d_' SeriesDescription '_lasttimept_' SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID
DCE_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' DCE_file_infix '.nii.gz']); 

for jj = 1:size(series_DCE_maps,1)
    SeriesInstanceUID = series_DCE_maps(jj).SeriesInstanceUID;
    SeriesDescription = strrep(strrep(strrep(strrep(strrep(series_DCE_maps(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    DCE_map_file_infix = sprintf(['%02d_' SeriesDescription '_' SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID 
    DCE_map_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' DCE_map_file_infix '.nii.gz']); 
    DCE_map_nii.hdr = DCE_nii.hdr; 
    save_nii(DCE_map_nii, [script03_prefix '00_radpath_raw/radpath_raw_' DCE_map_file_infix '.nii.gz'])
end

    
series_DSC = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%DSC%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);
series_DSC_maps = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%rBV%'' '...
        'or SeriesDescription like ''%rBF%'' '...
        'or SeriesDescription like ''%MTT%'' '...
        'or SeriesDescription like ''%Delay%'' '...
        'or SeriesDescription like ''%Leakage (K2)%'')'...
        'order by SeriesDate, SeriesInstanceUID ASC']);
SeriesInstanceUID = series_DSC.SeriesInstanceUID;
SeriesDescription = strrep(strrep(strrep(strrep(strrep(series_DSC.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
DSC_file_infix = sprintf(['%02d_' SeriesDescription '_lasttimept_' SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID
DSC_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' DSC_file_infix '.nii.gz']); 

for jj = 1:size(series_DSC_maps,1)
    SeriesInstanceUID = series_DSC_maps(jj).SeriesInstanceUID;
    SeriesDescription = strrep(strrep(strrep(strrep(strrep(series_DSC_maps(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    DSC_map_file_infix = sprintf(['%02d_' SeriesDescription '_' SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID 
    DSC_map_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' DSC_map_file_infix '.nii.gz']); 
    DSC_map_nii.hdr = DSC_nii.hdr; 
    save_nii(DSC_map_nii, [script03_prefix '00_radpath_raw/radpath_raw_' DSC_map_file_infix '.nii.gz'])
end

%% II. Register DCE/DSC lasttimept to Reference, then apply to all derived maps
tic
clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';

mksqlite('open', 'ctkDICOM.sql' ) ;
series_all = mksqlite(['select * from Series order by SeriesDate ASC']);
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);

series_T2FLAIR = mksqlite(['select * from Series where SeriesDescription like ''%FLAIR%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);
series_DCE = mksqlite(['select * from Series where SeriesDescription like ''%DCE%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);
series_DSC = mksqlite(['select * from Series where SeriesDescription like ''%DSC%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);

n_series_all = size(series_all,1);
n_studies_all = size(studies_all,1);
n_series_DCE = size(series_DCE,1);
n_series_DSC = size(series_DSC,1);
n_register = n_series_DCE + n_series_DSC; % Leave out the fixed image

%%%%%%%%%%%%%%%% 1. 
fid = fopen([script03_prefix 'Mass_Register.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_register
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');


zz=1;
for ii=4%1:n_studies_all
    fixed_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%FLAIR%'' ']);
    series_DCE = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%DCE%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);
    series_DCE_maps =     mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%K12%'' '...
        'or SeriesDescription like ''%K21%'' '...
        'or SeriesDescription like ''%Vp%'' '...
        'or SeriesDescription like ''%Ve (distribution volume)%'' '...
        'or SeriesDescription like ''%Wash-in%'' '...
        'or SeriesDescription like ''%Wash-out%'' '...
        'or SeriesDescription like ''%Time to peak%'' '...
        'or SeriesDescription like ''%Area under curve%'' '...
        'or SeriesDescription like ''%Peak enhancement%'')'...
        'order by SeriesDate, SeriesInstanceUID ASC']);
    
    n_series_DCE=size(series_DCE,1);
    n_series_DCE_maps = size(series_DCE_maps,1);
    
    Fixed_SeriesInstanceUID = fixed_series.SeriesInstanceUID;
    Moving_SeriesInstanceUID = series_DCE.SeriesInstanceUID;
    Fixed_SeriesDescription = strrep(strrep(strrep(strrep(strrep(fixed_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    Moving_SeriesDescription = strrep(strrep(strrep(strrep(strrep(series_DCE.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    Fixed_file_infix = sprintf(['%02d_' Fixed_SeriesDescription '_' Fixed_SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID
    Moving_file_infix = sprintf(['%02d_' Moving_SeriesDescription '_lasttimept_' Moving_SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID

    fprintf(fid,['\n\njob' num2str(zz) ':\n']);

    % 1. Perform rough down-sampled registration as beginning step to help
    % initialize later registration
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
        script03_prefix '08_RInput/Affine_Fixed_' Fixed_file_infix '.nii.gz 4 4 4 1\n']); % Down-sample fixed image
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix '.nii.gz '...
        script03_prefix '08_RInput/Affine_Moving_' Moving_file_infix '.nii.gz 4 4 4 1\n']); %Down-sample moving image
    fprintf(fid,['\t$(ANTSPATH)/antsAffineInitializer 3 '...
        script03_prefix '08_RInput/Affine_Fixed_' Fixed_file_infix '.nii.gz '...
        script03_prefix '08_RInput/Affine_Moving_' Moving_file_infix '.nii.gz '...
        script03_prefix '08_RInput/InitialAffine_' Moving_file_infix '.mat 15 0.1 0 10\n']); % Rough registration

    % 2. Create Laplacian (edge detection) version of images to assist in registration
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script03_prefix '08_RInput/Laplacian_' Fixed_file_infix '.nii.gz '...
        'Laplacian '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz \n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script03_prefix '08_RInput/Laplacian_' Moving_file_infix '.nii.gz '...
        'Laplacian '...
        script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix '.nii.gz 1.5 1\n']);

    % 3. Perform registration 
    fprintf(fid,['\t$(ANTSPATH)/antsRegistration -d 3 -u 1 -w [0.025, 0.975] '...
        '-o ' script03_prefix '09_ROutput/RegistrationOutput_' Moving_file_infix '_ '...
        '-r ' script03_prefix '08_RInput/InitialAffine_' Moving_file_infix '.mat -z 1 --float 0 '...
        '-m MI[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
              script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m MI[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
              script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
              '\n']);
%         '-m CC[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
%               script_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz,0.5,4] '...
%         '-m CC[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
%              script_prefix '08_RInput/Laplacian_' Moving_file_infix '.nii.gz,0.5,4] -c [50x10x0,1e-9,15] -t SyN[0.1,3,0] -f 4x2x1 -s 2x1x0'...
%                 '\n']);

    % 4. Apply registration 
    fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
        '-i ' script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix '.nii.gz '...
        '-o ' script03_prefix '10_Warped/Warped_' Moving_file_infix '.nii.gz '...
        '-r ' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
        '-n NearestNeighbor '... % Interpolation method for resampling
        '-t [' script03_prefix '09_ROutput/RegistrationOutput_' Moving_file_infix '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
        '\n']); 

%             '-t ' script_prefix '09_ROutput/RegistrationOutput_' Moving_file_infix '_1Warp.nii.gz --float 0'... % Deformation Field Transform
%             '\n']); 

    for jj = 1:n_series_DCE_maps
        Moving_SeriesInstanceUID_maps = series_DCE_maps(jj).SeriesInstanceUID;
        Moving_SeriesDescription_maps = strrep(strrep(strrep(strrep(strrep(series_DCE_maps(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
        Moving_file_infix_maps = sprintf(['%02d_' Moving_SeriesDescription_maps '_' Moving_SeriesInstanceUID_maps],ii); % patient #, Descrip, SeriesInstanceUID 
        fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
                '-i ' script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix_maps '.nii.gz '...
                '-o ' script03_prefix '10_Warped/Warped_' Moving_file_infix_maps '.nii.gz '... 
                '-r ' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
                '-n NearestNeighbor '... % Interpolation method for resampling
                '-t [' script03_prefix '09_ROutput/RegistrationOutput_' Moving_file_infix '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
                '\n']); 
    end
    
    zz=zz+1;
end

 
for ii=4%1:n_studies_all
    fixed_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%FLAIR%'' ']);
    series_DSC = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%DSC%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);
    series_DSC_maps = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%rBV%'' '...
        'or SeriesDescription like ''%rBF%'' '...
        'or SeriesDescription like ''%MTT%'' '...
        'or SeriesDescription like ''%Delay%'' '...
        'or SeriesDescription like ''%Leakage (K2)%'')'...
        'order by SeriesDate, SeriesInstanceUID ASC']);
    n_series_DSC=size(series_DSC,1);
    n_series_DSC_maps = size(series_DSC_maps,1);
    
    Fixed_SeriesInstanceUID = fixed_series.SeriesInstanceUID;
    Moving_SeriesInstanceUID = series_DSC.SeriesInstanceUID;
    Fixed_SeriesDescription = strrep(strrep(strrep(strrep(strrep(fixed_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    Moving_SeriesDescription = strrep(strrep(strrep(strrep(strrep(series_DSC.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    Fixed_file_infix = sprintf(['%02d_' Fixed_SeriesDescription '_' Fixed_SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID
    Moving_file_infix = sprintf(['%02d_' Moving_SeriesDescription '_lasttimept_' Moving_SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID

    fprintf(fid,['\n\njob' num2str(zz) ':\n']);

    % 1. Perform rough down-sampled registration as beginning step to help
    % initialize later registration
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
        script03_prefix '08_RInput/Affine_Fixed_' Fixed_file_infix '.nii.gz 4 4 4 1\n']); % Down-sample fixed image
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix '.nii.gz '...
        script03_prefix '08_RInput/Affine_Moving_' Moving_file_infix '.nii.gz 4 4 4 1\n']); %Down-sample moving image
    fprintf(fid,['\t$(ANTSPATH)/antsAffineInitializer 3 '...
        script03_prefix '08_RInput/Affine_Fixed_' Fixed_file_infix '.nii.gz '...
        script03_prefix '08_RInput/Affine_Moving_' Moving_file_infix '.nii.gz '...
        script03_prefix '08_RInput/InitialAffine_' Moving_file_infix '.mat 15 0.1 0 10\n']); % Rough registration

    % 2. Create Laplacian (edge detection) version of images to assist in registration
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script03_prefix '08_RInput/Laplacian_' Fixed_file_infix '.nii.gz '...
        'Laplacian '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz \n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script03_prefix '08_RInput/Laplacian_' Moving_file_infix '.nii.gz '...
        'Laplacian '...
        script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix '.nii.gz 1.5 1\n']);

    % 3. Perform registration 
    fprintf(fid,['\t$(ANTSPATH)/antsRegistration -d 3 -u 1 -w [0.025, 0.975] '...
        '-o ' script03_prefix '09_ROutput/RegistrationOutput_' Moving_file_infix '_ '...
        '-r ' script03_prefix '08_RInput/InitialAffine_' Moving_file_infix '.mat -z 1 --float 0 '...
        '-m MI[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
              script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m MI[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
              script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
              '\n']);
%         '-m CC[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
%               script_prefix '00_radpath_raw/' Moving_file_infix '.nii.gz,0.5,4] '...
%         '-m CC[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz,'...
%              script_prefix '08_RInput/Laplacian_' Moving_file_infix '.nii.gz,0.5,4] -c [50x10x0,1e-9,15] -t SyN[0.1,3,0] -f 4x2x1 -s 2x1x0'...
%                 '\n']);

    % 4. Apply registration 
    fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
        '-i ' script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix '.nii.gz '...
        '-o ' script03_prefix '10_Warped/Warped_' Moving_file_infix '.nii.gz '...
        '-r ' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
        '-n NearestNeighbor '... % Interpolation method for resampling
        '-t [' script03_prefix '09_ROutput/RegistrationOutput_' Moving_file_infix '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
        '\n']); 

%             '-t ' script_prefix '09_ROutput/RegistrationOutput_' Moving_file_infix '_1Warp.nii.gz --float 0'... % Deformation Field Transform
%             '\n']); 

    for jj = 1:n_series_DSC_maps
        Moving_SeriesInstanceUID_maps = series_DSC_maps(jj).SeriesInstanceUID;
        Moving_SeriesDescription_maps = strrep(strrep(strrep(strrep(strrep(series_DSC_maps(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');   
        Moving_file_infix_maps = sprintf(['%02d_' Moving_SeriesDescription_maps '_' Moving_SeriesInstanceUID_maps],ii); % patient #, Descrip, SeriesInstanceUID 
        fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
                '-i ' script03_prefix '00_radpath_raw/radpath_raw_' Moving_file_infix_maps '.nii.gz '...
                '-o ' script03_prefix '10_Warped/Warped_' Moving_file_infix_maps '.nii.gz '... 
                '-r ' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_SeriesInstanceUID '.nii.gz '...
                '-n NearestNeighbor '... % Interpolation method for resampling
                '-t [' script03_prefix '09_ROutput/RegistrationOutput_' Moving_file_infix '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
                '\n']); 
    end
    
    zz=zz+1;
end
   
fclose(fid);
system(['make -j 8 -f ' script03_prefix 'Mass_Register.makefile'])

time_s = toc
time_min = time_s/60
time_hrs = time_min/60 % 40 min

mksqlite('close' ) ;






%% Check results - make sure image volumes have no hiccups

clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';

mksqlite('open', 'ctkDICOM.sql' ) ;
series_all = mksqlite(['select * from Series order by SeriesDate ASC']);
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);

ii=4
fixed_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%FLAIR%'' ']);
temp_series = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%DCE%'' '...
        'or SeriesDescription like ''%K12%'' '...
        'or SeriesDescription like ''%K21%'' '...
        'or SeriesDescription like ''%Vp%'' '...
        'or SeriesDescription like ''%Ve (distribution volume)%'' '...
        'or SeriesDescription like ''%Wash-in%'' '...
        'or SeriesDescription like ''%Wash-out%'' '...
        'or SeriesDescription like ''%Time to peak%'' '...
        'or SeriesDescription like ''%Area under curve%'' '...
        'or SeriesDescription like ''%Peak enhancement%'''...
        'or SeriesDescription like ''%DSC%'' '...
        'or SeriesDescription like ''%rBV%'' '...
        'or SeriesDescription like ''%rBF%'' '...
        'or SeriesDescription like ''%MTT%'' '...
        'or SeriesDescription like ''%Delay%'' '...
        'or SeriesDescription like ''%Leakage (K2)%'') '...
        'order by SeriesDate, SeriesInstanceUID ASC']);  

jj=2;
SeriesInstanceUID = temp_series(jj).SeriesInstanceUID;
SeriesDescription = strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
file_infix = sprintf(['%02d_' SeriesDescription '_' SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID 


system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '07_HistMatch/Brain_HistMatch_' fixed_series.SeriesInstanceUID '.nii.gz '...
    '-o ' script03_prefix '10_Warped/Warped_' file_infix '.nii.gz'])

% for DCE/DSC lasttimept
file_infix = sprintf(['%02d_' SeriesDescription '_lasttimept_' SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID 
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap  '...
        '-g ' script03_prefix '00_radpath_raw/radpath_raw_' file_infix '.nii.gz']);
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap  '...
        '-g ' script03_prefix '10_Warped/Warped_' file_infix '.nii.gz']);
  
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '07_HistMatch/Brain_HistMatch_' fixed_series.SeriesInstanceUID '.nii.gz '...
    '-o ' script03_prefix '10_Warped/Warped_' file_infix '.nii.gz'])



system(['vglrun /opt/apps/SLICER/Slicer-4.3.1-linux-amd64/Slicer'])
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap'])

% View Histogram matched Brain-extracted volumes
system(['vglrun itksnap -g 00_radpath_raw/Brain_HistMatch_' SeriesInstanceUID '.nii.gz'])


