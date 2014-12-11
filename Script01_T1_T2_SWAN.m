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
mksqlite('open', 'ctkDICOM.sql' );
series_all = mksqlite(['select * from Series order by SeriesDate ASC']);
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);
n_series_all = size(series_all,1);
n_studies_all = size(studies_all,1);

fid = fopen([script01_prefix 'DICOM_to_NIfTI.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_series_all 
    fprintf(fid, [' job' num2str(ii)]);
end

zz=1;
for ii = 1:n_studies_all
    
    StudyInstanceUID = studies_all(ii).StudyInstanceUID;
    temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%T1%'' '... % Includes T1, T1post
        'or SeriesDescription like ''%T2%'' '... % Includes T2, T2star, T2FLAIR
        'or SeriesDescription like ''%SWAN%'') '... % Includes SWAN Image volume
        'order by SeriesNumber ASC']);
    numel(temp_series)
    for jj = 1:numel(temp_series) % will atuomatically pass over if no content in temp_series
        UID = temp_series(jj).SeriesInstanceUID;
        Descrip = strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
        ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); % patient #, Descrip, SeriesInstanceUID
        fprintf(fid,['\n\njob' num2str(zz) ':\n']);

        images   = mksqlite(['select * from Images where SeriesInstanceUID = ''' UID '''' ]);
        pieces = strsplit(images(1).Filename,'/'); % Only works in Matlab 2014a
        halves = strsplit(images(1).Filename,pieces(end)); % Use end fragment to figure out pathname
        pathname = halves{1,1};

        % Convert 1 series to Nifti - Command Source_Directory Output_Filename SeriesInstanceUID(to ensure only that series is used)
        fprintf(fid, ['\tDicomSeriesReadImageWrite2 ' pathname ' '...
            script01_prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz '...
            UID '\n']);
        zz=zz+1;
    end
end

fclose(fid);
system(['make -j 8 -f ' script01_prefix 'DICOM_to_NIfTI.makefile'])

mksqlite('close' ) ;






%% II. Perform Histogram Matching via these steps:
% 1. N4 Bias Correction
% 2. Brain Extraction
% 3. Intensity Truncation
% 4. Histogram Matching
%
tic
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
series_SWAN = mksqlite(['select * from Series where SeriesDescription like ''%SWAN%'' order by SeriesDate ASC']);%14, none for patient 7
series_T1 = mksqlite(['select * from Series where SeriesDescription like ''%T1%'' and (SeriesDescription not like ''%STEALTH%'' and SeriesDescription not like ''%+C%'') order by SeriesDate ASC']); %17
series_T1post = mksqlite(['select * from Series where SeriesDescription like ''%T1%'' and (SeriesDescription like ''%STEALTH%'' or SeriesDescription like ''%+C%'') order by SeriesDate ASC']); %17
series_T2 = mksqlite(['select * from Series where SeriesDescription like ''%T2%'' and SeriesDescription not like ''%*%'' and SeriesDescription not like ''%FLAIR%'' order by SeriesDate ASC']); %17
series_T2star = mksqlite(['select * from Series where SeriesDescription like ''%T2*%'' order by SeriesDate ASC']); %17
series_T2FLAIR = mksqlite(['select * from Series where SeriesDescription like ''%FLAIR%'' order by SeriesDate ASC']); %17

n_SWAN = size(series_SWAN,1);
n_T1 = size(series_T1,1);
n_T1post = size(series_T1post,1);
n_T2 = size(series_T2,1);
n_T2star = size(series_T2star,1);
n_T2FLAIR = size(series_T2FLAIR,1);


%% Histogram Match for the 1st set
%%%%%%%%%%%%%%%%%%%%%% 1.Repeat for T1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen([script01_prefix  'T1_Create_Mask.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_T1 
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:n_T1 
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    Create_Mask(fid, ii, series_T1, 5, 4);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'T1_Create_Mask.makefile'])



for ii=1:n_T1 
    UID = series_T1(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T1(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %

    % 2i. Extract Brain
    system(['\t$(ANTSPATH)/MultiplyImages 3 '...
		script01_prefix '04_Masks/Atropos_LC_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']); 
end


for ii=1:n_T1
    UID = series_T1(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T1(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    % 3a. Create cumulative histogram
    if ii==1 % For 1st image, initialize the XX_cume variable
        T1_cume = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        T1_cume.img = reshape(T1_cume.img,1,[]);
    else
        T1_temp = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        T1_cume.img = [T1_cume.img reshape(T1_temp.img,1,[])];
    end
    
    
end
imgsize_factors = factor(numel(T1_cume.img));
n3 = numel(imgsize_factors);
n1 = round(n3*1/3);
n2 = round(n3*2/3);
cume_dim1 = prod(imgsize_factors(1:n1));
cume_dim2 = prod(imgsize_factors(n1+1:n2));
cume_dim3 = prod(imgsize_factors(n2+1:n3));
T1_cume.hdr.dime.dim = [3 cume_dim1 cume_dim2 cume_dim3 1 1 1 1];
T1_cume.img = reshape(T1_cume.img,cume_dim1,cume_dim2,cume_dim3);
save_nii(T1_cume, 'T1_cume.nii.gz')
clear T1_cume T1_temp

% 3b. Histogram Match
fid = fopen([script01_prefix 'T1_HistMatch.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_T1 
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');
for ii=1:n_T1 
    UID = series_T1(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T1(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' ptno_Descrip_UID '.nii.gz '...
        'HistogramMatch '...
        script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
        script01_prefix 'T1_cume.nii.gz 255 64']);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'T1_HistMatch.makefile'])





%%
%%%%%%%%%%%%%%%%%%%%%% 2.Repeat for T1post %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen([script01_prefix 'T1post_Create_Mask.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_T1post
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:n_T1post
	fprintf(fid,['\n\njob' num2str(ii) ':\n']);
	Create_Mask(fid, ii, series_SWAN, 5, 4);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'T1post_Create_Mask.makefile'])



for ii=1:n_T1 
    UID = series_T1post(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T1post(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %

    % 2i. Extract Brain
    system(['\t$(ANTSPATH)/MultiplyImages 3 '...
		script01_prefix '04_Masks/Atropos_LC_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']); 
end
    
    
for ii=1:n_T1post
    UID = series_T1post(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T1post(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    % 3a. Create cumulative histogram
    if ii==1 % For 1st image, initialize the XX_cume variable
        T1post_cume = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        T1post_cume.img = reshape(T1post_cume.img,1,[]);
    else
        T1post_temp = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        T1post_cume.img = [T1post_cume.img reshape(T1post_temp.img,1,[])];
    end
end
imgsize_factors = factor(numel(T1post_cume.img));
n3 = numel(imgsize_factors);
n1 = round(n3*1/3);
n2 = round(n3*2/3);
cume_dim1 = prod(imgsize_factors(1:n1));
cume_dim2 = prod(imgsize_factors(n1+1:n2));
cume_dim3 = prod(imgsize_factors(n2+1:n3));
T1post_cume.hdr.dime.dim = [3 cume_dim1 cume_dim2 cume_dim3 1 1 1 1];
T1post_cume.img = reshape(T1post_cume.img,cume_dim1,cume_dim2,cume_dim3);
save_nii(T1post_cume, 'T1post_cume.nii.gz')
clear T1post_cume T1post_temp

% 3b. Histogram Match
fid = fopen([script01_prefix 'T1post_HistMatch.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_T1post
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');
for ii=1:n_T1post
    UID = series_T1post(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T1post(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' ptno_Descrip_UID '.nii.gz '...
        'HistogramMatch '...
        script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
        script01_prefix 'T1post_cume.nii.gz 255 64']);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'T1post_HistMatch.makefile'])


%%
%%%%%%%%%%%%%%%%%%%%%% 3.Repeat for T2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen([script01_prefix 'T2_Create_Mask.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_T2
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:n_T2
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    Create_Mask(fid, ii, series_T2, 5, 4);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'T2_Create_Mask.makefile'])
    

for ii=1:n_T1 
    UID = series_T2(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T2(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %

    % 2i. Extract Brain
    system(['\t$(ANTSPATH)/MultiplyImages 3 '...
		script01_prefix '04_Masks/Atropos_LC_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']); 
end


for ii=1:n_T2
    UID = series_T2(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T2(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    % 3a. Create cumulative histogram
    if ii==1 % For 1st image, initialize the XX_cume variable
        T2_cume = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        T2_cume.img = reshape(T2_cume.img,1,[]);
    else
        T2_temp = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        T2_cume.img = [T2_cume.img reshape(T2_temp.img,1,[])];
    end
end
imgsize_factors = factor(numel(T2_cume.img));
n3 = numel(imgsize_factors);
n1 = round(n3*1/3);
n2 = round(n3*2/3);
cume_dim1 = prod(imgsize_factors(1:n1));
cume_dim2 = prod(imgsize_factors(n1+1:n2));
cume_dim3 = prod(imgsize_factors(n2+1:n3));
T2_cume.hdr.dime.dim = [3 cume_dim1 cume_dim2 cume_dim3 1 1 1 1];
T2_cume.img = reshape(T2_cume.img,cume_dim1,cume_dim2,cume_dim3);
save_nii(T2_cume, 'T2_cume.nii.gz')
clear T2_cume T2_temp

% 3b. Histogram Match
fid = fopen([script01_prefix 'T2_HistMatch.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_T2
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');
for ii=1:n_T2
    UID = series_T2(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T2(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' ptno_Descrip_UID '.nii.gz '...
        'HistogramMatch '...
        script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
        script01_prefix 'T2_cume.nii.gz 255 64']);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'T2_HistMatch.makefile'])


%%
%%%%%%%%%%%%%%%%%%%%%% 4.Repeat for T2star %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen([script01_prefix 'T2star_Create_Mask.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_T2star
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:n_T2star 
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    Create_Mask(fid, ii, series_T2star, 5, 4);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'T2star_Create_Mask.makefile'])
    




for ii=1:n_T1 
    UID = series_T2star(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T2star(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %

    % 2i. Extract Brain
    system(['\t$(ANTSPATH)/MultiplyImages 3 '...
		script01_prefix '04_Masks/Atropos_LC_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']); 
end

    

for ii=1:n_T2star
    UID = series_T2star(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T2star(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star'); 
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    % 3a. Create cumulative histogram
    if ii==1 % For 1st image, initialize the XX_cume variable
        T2star_cume = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        T2star_cume.img = reshape(T2star_cume.img,1,[]);
    else
        T2star_temp = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        T2star_cume.img = [T2star_cume.img reshape(T2star_temp.img,1,[])];
    end
end
imgsize_factors = factor(numel(T2star_cume.img));
n3 = numel(imgsize_factors);
n1 = round(n3*1/3);
n2 = round(n3*2/3);
cume_dim1 = prod(imgsize_factors(1:n1));
cume_dim2 = prod(imgsize_factors(n1+1:n2));
cume_dim3 = prod(imgsize_factors(n2+1:n3));
T2star_cume.hdr.dime.dim = [3 cume_dim1 cume_dim2 cume_dim3 1 1 1 1];
T2star_cume.img = reshape(T2star_cume.img,cume_dim1,cume_dim2,cume_dim3);
save_nii(T2star_cume, 'T2star_cume.nii.gz')
clear T2star_cume T2star_temp

% 3b. Histogram Match
fid = fopen([script01_prefix 'T2star_HistMatch.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_T2star
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');
for ii=1:n_T2star
    UID = series_T2star(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T2star(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star'); 
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' ptno_Descrip_UID '.nii.gz '...
        'HistogramMatch '...
        script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
        script01_prefix 'T2star_cume.nii.gz 255 64']);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'T2star_HistMatch.makefile'])


%%
%%%%%%%%%%%%%%%%%%%%%% 5.Repeat for T2FLAIR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen([script01_prefix 'T2FLAIR_Create_Mask.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_T2FLAIR 
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:n_T2FLAIR 
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    Create_Mask(fid, ii, series_T2FLAIR, 5, 4);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'T2FLAIR_Create_Mask.makefile'])


for ii=1:n_T1 
    UID = series_T2FLAIR(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T2FLAIR(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %

    % 2i. Extract Brain
    system(['\t$(ANTSPATH)/MultiplyImages 3 '...
		script01_prefix '04_Masks/Atropos_LC_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']); 
end

for ii=1:n_T2FLAIR
    UID = series_T2FLAIR(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T2FLAIR(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');  
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    % 3a. Create cumulative histogram
    if ii==1 % For 1st image, initialize the XX_cume variable
        T2FLAIR_cume = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        T2FLAIR_cume.img = reshape(T2FLAIR_cume.img,1,[]);
    else
        T2FLAIR_temp = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        T2FLAIR_cume.img = [T2FLAIR_cume.img reshape(T2FLAIR_temp.img,1,[])];
    end
end
imgsize_factors = factor(numel(T2FLAIR_cume.img));
n3 = numel(imgsize_factors);
n1 = round(n3*1/3);
n2 = round(n3*2/3);
cume_dim1 = prod(imgsize_factors(1:n1));
cume_dim2 = prod(imgsize_factors(n1+1:n2));
cume_dim3 = prod(imgsize_factors(n2+1:n3));
T2FLAIR_cume.hdr.dime.dim = [3 cume_dim1 cume_dim2 cume_dim3 1 1 1 1];
T2FLAIR_cume.img = reshape(T2FLAIR_cume.img,cume_dim1,cume_dim2,cume_dim3);
save_nii(T2FLAIR_cume, 'T2FLAIR_cume.nii.gz')
clear T2FLAIR_cume T2FLAIR_temp

% 3b. Histogram Match
fid = fopen([script01_prefix 'T2FLAIR_HistMatch.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_T2FLAIR
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');
for ii=1:n_T2FLAIR
    UID = series_T2FLAIR(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_T2FLAIR(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');    
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' ptno_Descrip_UID '.nii.gz '...
        'HistogramMatch '...
        script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
        script01_prefix 'T2FLAIR_cume.nii.gz 255 64']);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'T2FLAIR_HistMatch.makefile'])

        
time_s = toc
time_min = time_s/60
time_hrs = time_min/60 % 33 hours, next time use 5 cores?

mksqlite('close' ) ;

%%
%%%%%%%%%%%%%%%% 6. Perform Histogram Match for SWAN   %%%%%%%%%%%%%%
fid = fopen([script01_prefix 'SWAN_Create_Mask.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_SWAN 
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

zz=1; % Counter for every single job - combine into 1 long file?
for ii=4%1:n_SWAN 
    fprintf(fid,['\n\njob' num2str(zz) ':\n']);
    Create_Mask(fid, ii, series_SWAN, 5, 4);
    zz = zz+1;
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'SWAN_Create_Mask.makefile'])

% Check mask - how well it matches head volume
system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '00_radpath_raw/radpath_raw_' UID '.nii.gz '])


% Use mask to extract brain from images
for ii=1:n_T1 

    UID = series_SWAN(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_SWAN(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %

    % 2i. Extract Brain
    system(['\t$(ANTSPATH)/MultiplyImages 3 '...
		script01_prefix '04_Masks/Atropos_LC_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']); 
end

% Create cumulative histogram
for ii=1:n_SWAN
    UID = series_SWAN(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_SWAN(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    
    % 3a. Create cumulative histogram
    if ii==1 % For 1st SWAN image, initialize the SWAN_cume variable
        SWAN_cume = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' UID '.nii.gz']);
        SWAN_cume.img = reshape(SWAN_cume.img,1,[]);
    else
        SWAN_temp = load_nii([script01_prefix '06_Truncated/Brain_Truncated0_' UID '.nii.gz']);
        SWAN_cume.img = [SWAN_cume.img reshape(SWAN_temp.img,1,[])];
    end
end
imgsize_factors = factor(numel(SWAN_cume.img));
n3 = numel(imgsize_factors);
n1 = round(n3*1/3);
n2 = round(n3*2/3);
cume_dim1 = prod(imgsize_factors(1:n1));
cume_dim2 = prod(imgsize_factors(n1+1:n2));
cume_dim3 = prod(imgsize_factors(n2+1:n3));
SWAN_cume.hdr.dime.dim = [3 cume_dim1 cume_dim2 cume_dim3 1 1 1 1];
SWAN_cume.img = reshape(SWAN_cume.img,cume_dim1,cume_dim2,cume_dim3);
save_nii(SWAN_cume, [script01_prefix 'SWAN_cume.nii.gz'])
clear SWAN_cume SWAN_temp

% 3b. Histogram Match
fid = fopen([script01_prefix 'SWAN_HistMatch.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_SWAN 
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');
for ii=1:n_SWAN 
    UID = series_SWAN(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_SWAN(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' ptno_Descrip_UID '.nii.gz '...
        'HistogramMatch '...
        script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
        script01_prefix 'SWAN_cume.nii.gz 255 64']);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'SWAN_HistMatch.makefile'])










%% III. Register each image volume to reference for the same patient
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
series_SWAN = mksqlite(['select * from Series where SeriesDescription like ''%SWAN%'' order by SeriesDate ASC']);%14, none for patient 6 or 7
series_T1 = mksqlite(['select * from Series where SeriesDescription like ''%Ax T1%'' and SeriesDescription not like ''%+C%'' order by SeriesDate ASC']); %15
series_T1post = mksqlite(['select * from Series where SeriesDescription like ''%T1%'' and SeriesDescription like ''%STEALTH%'' or SeriesDescription like ''%+C%'' order by SeriesDate ASC']); %15
series_T2 = mksqlite(['select * from Series where SeriesDescription like ''%T2%'' and SeriesDescription not like ''%*%'' and SeriesDescription not like ''%FLAIR%'' order by SeriesDate ASC']); %15
series_T2star = mksqlite(['select * from Series where SeriesDescription like ''%T2*%'' order by SeriesDate ASC']); %15
series_T2FLAIR = mksqlite(['select * from Series where SeriesDescription like ''%FLAIR%'' order by SeriesDate ASC']); %15

n_series_all = size(series_all,1);
n_studies_all = size(studies_all,1);
n_SWAN = size(series_SWAN,1);
n_T1 = size(series_T1,1);
n_T1post = size(series_T1post,1);
n_T2 = size(series_T2,1);
n_T2star = size(series_T2star,1);
n_T2FLAIR = size(series_T2FLAIR,1);

n_register = n_SWAN + n_T1 + n_T1post + n_T2 + n_T2star; % Leave out the fixed image

%%%%%%%%%%%%%%%% 1. 
fid = fopen([script01_prefix 'Mass_Register.makefile'],'w');
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
        'and (SeriesDescription like ''%SWAN%'' '... % Includes SWAN Image volume
        'or SeriesDescription like ''%T1%'' '... % Includes T1, T1post
        'or SeriesDescription like ''%T2%'' '... % Includes T2, T2star, T2FLAIR
        'or SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
        'or SeriesDescription like ''%Fractional Aniso%'' '... % Includes FA from DTI
        'or SeriesDescription like ''%Average DC%'') '...  % Includes ADC from DTI
        'and SeriesDescription not like ''%FLAIR%'' '... % Exclude the fixed image volume
        'order by SeriesNumber ASC']);
    n_moving_series=size(moving_series,1);
    
    for jj=1:n_moving_series
        Fixed_SeriesInstanceUID = fixed_series.SeriesInstanceUID;
        Moving_SeriesInstanceUID = moving_series(jj).SeriesInstanceUID;
        Fixed_SeriesDescription = strrep(strrep(strrep(strrep(strrep(fixed_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
        Moving_SeriesDescription = strrep(strrep(strrep(strrep(strrep(moving_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');        
        Fixed_file_infix = sprintf(['%02d_' Fixed_SeriesDescription '_' Fixed_SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID
        Moving_file_infix = sprintf(['%02d_' Moving_SeriesDescription '_' Moving_SeriesInstanceUID],ii); % patient #, Descrip, SeriesInstanceUID
        
        fprintf(fid,['\n\njob' num2str(zz) ':\n']);

        % 1. Perform rough down-sampled registration as beginning step to help
        % initialize later registration
        fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
            script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_file_infix '.nii.gz '...
            script01_prefix '08_RInput/Affine_Fixed_' Fixed_file_infix '.nii.gz 4 4 4 1\n']); % Down-sample fixed image
        fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
            script01_prefix '07_HistMatch/Brain_HistMatch_' Moving_file_infix '.nii.gz '...
            script01_prefix '08_RInput/Affine_Moving_' Moving_file_infix '.nii.gz 4 4 4 1\n']); %Down-sample moving image
        fprintf(fid,['\t$(ANTSPATH)/antsAffineInitializer 3 '...
            script01_prefix '08_RInput/Affine_Fixed_' Fixed_file_infix '.nii.gz '...
            script01_prefix '08_RInput/Affine_Moving_' Moving_file_infix '.nii.gz '...
            script01_prefix '08_RInput/InitialAffine_' Moving_file_infix '.mat 15 0.1 0 10\n']); % Rough registration

        % 2. Create Laplacian (edge detection) version of images to assist in registration
        fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
            script01_prefix '08_RInput/Laplacian_' Fixed_file_infix '.nii.gz '...
            'Laplacian '...
            script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_file_infix '.nii.gz 1.5 1\n']);
        fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
            script01_prefix '08_RInput/Laplacian_' Moving_file_infix '.nii.gz '...
            'Laplacian '...
            script01_prefix '07_HistMatch/Brain_HistMatch_' Moving_file_infix '.nii.gz 1.5 1\n']);

        % 3. Perform registration 
        fprintf(fid,['\t$(ANTSPATH)/antsRegistration -d 3 -u 1 -w [0.025, 0.975] '...
            '-o ' script01_prefix '09_ROutput/RegistrationOutput_' Moving_file_infix '_ '...
            '-r ' script01_prefix '08_RInput/InitialAffine_' Moving_file_infix '.mat -z 1 --float 0 '...
            '-m MI[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_file_infix '.nii.gz,'...
                  script01_prefix '07_HistMatch/Brain_HistMatch_' Moving_file_infix '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
            '-m MI[' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_file_infix '.nii.gz,'...
                  script01_prefix '07_HistMatch/Brain_HistMatch_' Moving_file_infix '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
                  '> StandInForNameRegOutputMetric.txt\n']);
%         '-m CC[07_HistMatch/Brain_HistMatch_' Fixed_file_infix '.nii.gz,'...
%               '07_HistMatch/Brain_HistMatch_' Moving_file_infix '.nii.gz,0.5,4] '...
%         '-m CC[08_RInput/Laplacian_' Fixed_file_infix '.nii.gz,'...
%              '08_RInput/Laplacian_' Moving_file_infix '.nii.gz,0.5,4] -c [50x10x0,1e-9,15] -t SyN[0.1,3,0] -f 4x2x1 -s 2x1x0'...
%                 '\n']);
        
        % 4. Apply registration 
        fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
            '-i ' script01_prefix '07_HistMatch/Brain_HistMatch_' Moving_file_infix '.nii.gz '...
            '-o ' script01_prefix '10_Warped/Warped_' Moving_file_infix '.nii.gz '...
            '-r ' script01_prefix '07_HistMatch/Brain_HistMatch_' Fixed_file_infix '.nii.gz '...
            '-n NearestNeighbor '... % Interpolation method for resampling
            '-t [' script01_prefix '09_ROutput/RegistrationOutput_' Moving_file_infix '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
            '\n']); 
            
%             '-t 09_ROutput/RegistrationOutput_' Moving_file_infix '_1Warp.nii.gz --float 0'... % Deformation Field Transform
%             '\n']); 

        zz=zz+1;
    end
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'Mass_Register.makefile'])

time_s = toc
time_min = time_s/60
time_hrs = time_min/60 % 33 hours, next time use 5 cores?

mksqlite('close' ) ;






%%


% Check results - make sure image volumes have no hiccups

ii=2; % Patient #
SWAN_UID = series_SWAN(ii).SeriesInstanceUID; % good, except #3
SWAN_Descrip = strrep(strrep(strrep(strrep(strrep(series_SWAN(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
SWAN_ptno_Descrip_UID = sprintf(['%02d_' SWAN_Descrip '_' SWAN_UID],ii); %
    
T1_UID = series_T1(ii).SeriesInstanceUID; % good
T1_Descrip = strrep(strrep(strrep(strrep(strrep(series_T1(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
T1_ptno_Descrip_UID = sprintf(['%02d_' T1_Descrip '_' T1_UID],ii); %

T1post_UID = series_T1post(ii).SeriesInstanceUID; % small - misses top vessel
T1post_Descrip = strrep(strrep(strrep(strrep(strrep(series_T1post(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
T1post_ptno_Descrip_UID = sprintf(['%02d_' T1post_Descrip '_' T1post_UID],ii); %

T2_UID = series_T2(ii).SeriesInstanceUID; % small - misses borders
T2_Descrip = strrep(strrep(strrep(strrep(strrep(series_T2(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
T2_ptno_Descrip_UID = sprintf(['%02d_' T2_Descrip '_' T2_UID],ii); %

T2star_UID = series_T2star(ii).SeriesInstanceUID; % small -misses borders
T2star_Descrip = strrep(strrep(strrep(strrep(strrep(series_T2star(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');   
T2star_ptno_Descrip_UID = sprintf(['%02d_' T2star_Descrip '_' T2star_UID],ii); %

T2FLAIR_UID = series_T2FLAIR(ii).SeriesInstanceUID; % good
T2FLAIR_Descrip = strrep(strrep(strrep(strrep(strrep(series_T2FLAIR(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star'); 
T2FLAIR_ptno_Descrip_UID = sprintf(['%02d_' T2FLAIR_Descrip '_' T2FLAIR_UID],ii); %

% 
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '00_radpath_raw/radpath_raw_' T2_UID '.nii.gz '])
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '00_radpath_raw/radpath_raw_' T2star_UID '.nii.gz '])
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '00_radpath_raw/radpath_raw_' T2FLAIR_UID '.nii.gz '])






display('SWAN')
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '01_N4/Brain_N4_' SWAN_UID '.nii.gz '...
    '-s ' script01_prefix '04_Masks/Mask_largestcomponent_' SWAN_UID '.nii.gz'])

display('T1')
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '01_N4/Brain_N4_' T1_UID '.nii.gz '...
    '-s ' script01_prefix '04_Masks/Mask_largestcomponent_' T1_UID '.nii.gz'])

display('T1post')
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '01_N4/Brain_N4_' T1post_UID '.nii.gz '...
    '-s ' script01_prefix '04_Masks/Mask_largestcomponent_' T1post_UID '.nii.gz'])

display('T2')
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '01_N4/Brain_N4_' T2_UID '.nii.gz '...
    '-s ' script01_prefix '04_Masks/Mask_largestcomponent_' T2_UID '.nii.gz'])

display('T2star')
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '01_N4/Brain_N4_' T2star_UID '.nii.gz '...
    '-s ' script01_prefix '04_Masks/Mask_largestcomponent_' T2star_UID '.nii.gz'])

display('T2FLAIR')
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '01_N4/Brain_N4_' T2FLAIR_UID '.nii.gz '...
    '-s ' script01_prefix '04_Masks/Mask_largestcomponent_' T2FLAIR_UID '.nii.gz'])


system(['vglrun /opt/apps/SLICER/Slicer-4.3.1-linux-amd64/Slicer'])

%% View Histogram matched Brain-extracted volumes
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
n_series_all = size(series_all,1);
n_studies_all = size(studies_all,1);

ii=4;
temp_series = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
    'and (SeriesDescription like ''%T1%'' '...
    'or SeriesDescription like ''%T2%'' '...
    'or SeriesDescription like ''%SWAN%'') '...
    'order by SeriesDate ASC']); 

jj=1;
UID = temp_series(jj).SeriesInstanceUID;
Descrip = strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/'...
    'ImageMath 3 ' script01_prefix 'Diff_' Descrip '.nii.gz  - '...
    script01_prefix '00_radpath_raw/radpath_raw_' UID '.nii.gz '...
    script01_prefix '07_HistMatch/Brain_HistMatch_' UID '.nii.gz '])

system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '00_radpath_raw/radpath_raw_' UID '.nii.gz '])
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '07_HistMatch/Brain_HistMatch_' UID '.nii.gz '])
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix 'Diff_' Descrip '.nii.gz '])


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

    
% Get histogram of cumulative one as well






% For pix of examples of metadata stored in SQL databse
studies_1 = studies_all(1);
series_1 = series_all(1);
images_all = mksqlite(['select * from Images where SeriesInstanceUID = ''' UID ''' ']); 
images_1 = images_all(1);




