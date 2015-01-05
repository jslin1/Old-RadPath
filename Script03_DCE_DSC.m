


clear all
close all
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
series_all = mksqlite(['select * from Series order by SeriesDate, SeriesInstanceUID ASC']);
series_DCE_DSC = mksqlite(['select * from Series '...
    'where SeriesDescription like ''%DCE%'' '...
    'or SeriesDescription like ''%DSC%'' '...
    'order by SeriesDate, SeriesInstanceUID ASC']);

    series_T2 = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%T2%'' '... % Includes T2, T2star, T2FLAIR
        'and SeriesDescription not like ''%*%'' '... % Exclude T2star
        'and SeriesDescription not like ''%FLAIR%'' '... % Exclude T2 FLAIR
        'order by SeriesDate ASC']);

    series_DCE = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%DCE%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);

    series_DCE_maps =     mksqlite(['select * from Series where '...
        'SeriesDescription like ''%K12%'' '...
        'or SeriesDescription like ''%K21%'' '...
        'or SeriesDescription like ''%Vp%'' '...
        'or SeriesDescription like ''%Ve (distribution volume)%'' '...
        'or SeriesDescription like ''%Wash-in%'' '...
        'or SeriesDescription like ''%Wash-out%'' '...
        'or SeriesDescription like ''%Time to peak%'' '...
        'or SeriesDescription like ''%Area under curve%'' '...
        'or SeriesDescription like ''%Peak enhancement%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);
    
    series_DSC = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%DSC%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);

    series_DSC_maps = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%rBV%'' '...
        'or SeriesDescription like ''%rBF%'' '...
        'or SeriesDescription like ''%MTT%'' '...
        'or SeriesDescription like ''%Delay%'' '...
        'or SeriesDescription like ''%Leakage (K2)%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);
    
    series_maps = mksqlite(['select * from Series where '...
            'SeriesDescription like ''%rBV%'' '...
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
            'or SeriesDescription like ''%Peak enhancement%'' '...
            'order by SeriesDate, SeriesInstanceUID ASC']);
        
series_DCE_DSC_maps = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%DCE%'' '...
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
        'or SeriesDescription like ''%Peak enhancement%'' '...
        'order by SeriesDate, SeriesInstanceUID ASC']);

prefix = script03_prefix;
series_T2_register = series_T2([2:8 10:size(series_T2,1)]); % T2+C

%% I. Re-number DICOMs for proper loading
dicom_series = series_DCE_DSC_maps;
for jj=1:size(dicom_series,1) % for NICE/NNL, Olea DROs
    jj
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(dicom_series(jj));
    images   = mksqlite(['select * from Images where SeriesInstanceUID = ''' UID ''' order by SOPInstanceUID ASC']);
    pieces = strsplit(images(1).Filename,'/'); % Only works in Matlab 2014a
    halves = strsplit(images(1).Filename,pieces(end)); % Use end fragment to figure out pathname
    pathname1 = halves{1,1}; % where it's stored
    pathname2 = [script03_prefix '00_radpath_raw/DICOMs_Re-numbered/' ptno '/' ptno '_' Descrip '_DICOMs' ]; % where it's going
    if ~isdir([pathname2])
        mkdir([pathname2]) % Create new directory (if it doesn't already exist)
    end
    
    % Generate list of DICOM filenames 
    if size(dir([pathname1 '*.dcm']),1)>0
        tmp_pathname = dir([pathname1 '*.dcm']);
    elseif size(dir([pathname1 '*IM*']),1)>0
        tmp_pathname = dir([pathname1 '*IM*']);
    end
    file_names = strvcat(tmp_pathname(1:size(tmp_pathname,1)).name); % Generates file names for every file in that directory
    nimages = size(file_names,1); % # of files in that directory
    
    % Copy DICOMs and re-name based on Instance Number 
    if ~isempty(strfind(Descrip,'DCE')) || ~isempty(strfind(Descrip,'DSC')) % if it's DCE or DSC, extract last timept, and re-number
        for ii = nimages-200:nimages; 
            tmp2=dicominfo(images(ii).Filename,'dictionary','gems-dicom-dict.txt');
            if tmp2.NumberOfTemporalPositions==tmp2.TemporalPositionIdentifier % if last timept (last timept == timept of current image)
                status = copyfile(images(ii).Filename, [pathname2 '/' sprintf('%05d',tmp2.InstanceNumber)],'f'); 
            end
        end
    else % if it's a map, then copy file, and re-name based on Instance Number
        for ii = 1:nimages; 
            tmp2=dicominfo(images(ii).Filename,'dictionary','gems-dicom-dict.txt');
            status = copyfile(images(ii).Filename, [pathname2 '/' sprintf('%05d',tmp2.InstanceNumber)],'f'); % finish this command
        end
    end
end


%% II. DICOM --> Nifti
n_cores = 9;
makefile = 'DICOM_to_NIfTI.makefile';
fid = fopen([prefix makefile],'w');
fprintf(fid, 'all:');
for ii=1:size(series_DCE_DSC_maps,1)
    fprintf(fid, [' job' num2str(ii)]);
end

dicom_series = series_DCE_DSC_maps;
for jj = 1:size(series_DCE_DSC_maps,1)
    fprintf(fid,['\n\njob' num2str(jj) ':\n']);
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_DCE_DSC_maps(jj));     
    pathname = [prefix '00_radpath_raw/DICOMs_Re-numbered/' ptno '/' ptno '_' Descrip '_DICOMs' ];    
    
    fprintf(fid, ['\tDicomSeriesReadImageWrite2 '...
    pathname ' '...
    prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz '...
    UID '\n']);
end

fclose(fid);
system(['make -j ' num2str(n_cores) ' -f ' prefix makefile])


%% III. Fix headers of DCE/DSC maps - spatial location info = applying registration works
for ii=1:size(series_DCE,1)
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_DCE(ii));
    DCE_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz']); 
    [field study]= find(cell2mat(cellfun(@(x) isequal(x,series_DCE(ii).StudyInstanceUID),struct2cell(series_DCE_maps),'UniformOutput',false)));
    for jj = study'
        [map_ptno map_Descrip map_UID map_ptno_Descrip_UID] = Generate_Label(series_DCE_maps(jj)); 
        map_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' map_ptno_Descrip_UID '.nii.gz']); 
        map_nii.hdr.dime.pixdim = DCE_nii.hdr.dime.pixdim;
        map_nii.hdr.hist = DCE_nii.hdr.hist;
        save_nii(map_nii, [map_nii.fileprefix '.nii.gz'])
    end
end

for ii=1:size(series_DSC,1)
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_DSC(ii));
    DSC_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz']); 
    [field study]= find(cell2mat(cellfun(@(x) isequal(x,series_DSC(ii).StudyInstanceUID),struct2cell(series_DSC_maps),'UniformOutput',false)));
    for jj = study'
        [map_ptno map_Descrip map_UID map_ptno_Descrip_UID] = Generate_Label(series_DSC_maps(jj)); 
        map_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' map_ptno_Descrip_UID '.nii.gz']); 
        map_nii.hdr.dime.pixdim = DSC_nii.hdr.dime.pixdim;
        map_nii.hdr.hist = DSC_nii.hdr.hist;
        save_nii(map_nii, [map_nii.fileprefix '.nii.gz'])
    end
end


%% IV. Register DCE/DSC lasttimept to Reference, then apply to all derived maps
tic
n_cores = 8;
fixed_series = series_T2_register;
Register3(prefix, 'DCE_register.makefile', n_cores, fixed_series, series_DCE, series_DCE_maps)
Register3(prefix, 'DSC_register.makefile', n_cores, fixed_series, series_DSC, series_DSC_maps)
time_s = toc
time_min = time_s/60
time_hrs = time_min/60 % 2.66 hrs


%% V. Create Masks
MD_parameter = 5;
MC_parameter = 4;
Create_Mask(prefix, 'DCE_Create_Mask.makefile', n_cores, series_DCE, MD_parameter, MC_parameter)
Create_Mask(prefix, 'DSC_Create_Mask.makefile', n_cores, series_DSC, MD_parameter, MC_parameter)


%% VI. Extract/Truncate
Extract_Truncate(prefix, series_DCE, series_DCE_maps)
Extract_Truncate(prefix, series_DSC, series_DSC_maps)


%% VII. Generate Cumulative Histogram - Further stratify by AIF type?
[rBV_MCA_R1_cume] = Generate_Cume(series_rBV_MCA_R1);
save_nii(rBV_MCA_R1_cume,[prefix 'rBV_MCA_R1_cume.nii.gz'])
clear rBV_MCA_R1_cume

%% VIII. Histogram Matching
n_intensity_levels = 255;
n_landmarks = 64;
HistMatch(prefix, 'T1_HistMatch.makefile',      'T1_cume.nii.gz',      n_cores, series_T1,      n_intensity_levels, n_landmarks)


%% Check results - make sure image volumes have no hiccups
jj=1;
[f_ptno f_Descrip f_UID f_ptno_Descrip_UID] = Generate_Label(fixed_series(jj));
[ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_DCE(jj));
system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '01_N4/N4_' f_ptno_Descrip_UID '.nii.gz '])

system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz '])


system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '01_N4/N4_' f_ptno_Descrip_UID '.nii.gz '...
    '-o ' script03_prefix '04_Warped/Warped_' ptno_Descrip_UID '.nii.gz'])

system(['vglrun /opt/apps/SLICER/Slicer-4.4.0-linux-amd64/Slicer'])




