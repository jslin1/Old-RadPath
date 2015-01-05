function Extract_Truncate(prefix, mask_series, extract_series)

for ii=1:size(mask_series,1)

[m_ptno m_Descrip m_UID m_ptno_Descrip_UID] = Generate_Label(mask_series(ii));
[e_field e_study]= find(cell2mat(cellfun(@(x) isequal(x,mask_series(ii).StudyInstanceUID),struct2cell(extract_series),'UniformOutput',false)));

for jj = e_study'
    [e_ptno e_Descrip e_UID e_ptno_Descrip_UID] = Generate_Label(extract_series(jj));
    % 1. Extract Brain
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/MultiplyImages 3 '...
		prefix '07_Masks/FINALMASK_' m_ptno_Descrip_UID '.nii.gz '...
		prefix '04_Warped/Warped_' e_ptno_Descrip_UID '.nii.gz '...
		prefix '08_Brain_Only/Brain_' e_ptno_Descrip_UID '.nii.gz\n']);

    % 2. Intensity Truncation
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/ImageMath 3 '...
		prefix '09_Truncated/Brain_Truncated0_' e_ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		prefix '08_Brain_Only/Brain_' e_ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']); 
end


end






































%%
%%%%%%%%%%%%%%%%%%%%%% 2.Repeat for T1post %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:n_T1post 
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_T1post(ii));

    % 2i. Extract Brain
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/MultiplyImages 3 '...
		script01_prefix '04_Masks/FINALMASK_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']); 
end
    
    
[T1post_cume] = Generate_Cume(series_T1post);
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
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_T1post(ii));
    
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
    

for ii=1:n_T2 
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_T2(ii));

    % 2i. Extract Brain
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/MultiplyImages 3 '...
		script01_prefix '04_Masks/FINALMASK_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']);
end


[T2_cume] = Generate_Cume(series_T2);
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
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_T2(ii));
    
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
    
for ii=1:n_T2star 
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_T2star(ii));

    % 2i. Extract Brain
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/MultiplyImages 3 '...
		script01_prefix '04_Masks/FINALMASK_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']);
end

 

[T2star_cume] = Generate_Cume(series_T2star);
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
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_T2star(ii));
    
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


for ii=1:n_T2FLAIR 
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_T2FLAIR(ii));

    % 2i. Extract Brain
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/MultiplyImages 3 '...
		script01_prefix '04_Masks/FINALMASK_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']);
end

[T2FLAIR_cume] = Generate_Cume(series_T2FLAIR);
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
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_T2FLAIR(ii));
    
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



% Check mask - how well it matches head volume    
UID = series_SWAN(ii).SeriesInstanceUID;
Descrip = strrep(strrep(strrep(strrep(strrep(series_SWAN(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); % patient #, Descrip, SeriesInstanceUID
system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz '...
    '-s ' script01_prefix '04_Masks/Atropos_' ptno_Descrip_UID '.nii.gz '])


% Use mask to extract brain from images
for ii=1:n_SWAN 

    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_SWAN(ii));

    % 2i. Extract Brain
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/MultiplyImages 3 '...
		script01_prefix '04_Masks/FINALMASK_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz\n']);

    % 3. Intensity Truncation
    system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/ImageMath 3 '...
		script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
		'TruncateImageIntensity '...
		script01_prefix '05_Brain_Only/Brain_' ptno_Descrip_UID '.nii.gz 0.01 0.999 256\n']);
end

% Create cumulative histogram
[SWAN_cume] = Generate_Cume(series_SWAN);
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
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_SWAN(ii));

    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        script01_prefix '07_HistMatch/Brain_HistMatch_' ptno_Descrip_UID '.nii.gz '...
        'HistogramMatch '...
        script01_prefix '06_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
        script01_prefix 'SWAN_cume.nii.gz 255 64']);
end
fclose(fid);
system(['make -j 8 -f ' script01_prefix 'SWAN_HistMatch.makefile'])










%%


% Check results - make sure image volumes have no hiccups
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


[ptno Descrip UID ptno_Descrip_UID] = Generate_Label(temp_series(jj));
       


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

ii=4;
temp_series = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
    'and (SeriesDescription like ''%T1%'' '...
    'or SeriesDescription like ''%T2%'' '...
    'or SeriesDescription like ''%SWAN%'') '...
    'order by SeriesDate ASC']); 

jj=1;
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(temp_series(jj));
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




