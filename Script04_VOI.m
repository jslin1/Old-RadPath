%% I. Load stuff
cd('/mnt/data/scratch/igilab/jslin1/RadPath')
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
series_T1 = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%Ax T1%'' '... % Includes T1, T1post
    'and SeriesDescription not like ''%+C%'' '... % Excludes T1post
    'order by SeriesDate ASC']); 
series_T1post = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%T1%'' '... % Includes T1, T1post
    'and (SeriesDescription like ''%STEALTH%'' '... %Includes only post
    'or SeriesDescription like ''%+C%'') '... % Includes only post
    'order by SeriesNumber ASC']);
series_T2 = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%T2%'' '... % Includes T2, T2star, T2 FLAIR
    'and SeriesDescription not like ''%*%'' '... % Excludes T2star
    'and SeriesDescription not like ''%FLAIR%'' '... % Excludes T2 FLAIR
    'order by SeriesDate ASC']);
series_T2star = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%T2*%'' '... % Includes T2star
    'order by SeriesDate ASC']);
series_T2FLAIR = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%FLAIR%'' '... % Includes T2FLAIR
    'order by SeriesDate ASC']);
series_SWAN = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%SWAN%'' '... % Includes SWAN
    'order by SeriesDate ASC']);

%% DWI/DTI
series_ADC = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC
        'and SeriesDescription not like ''%Exponential%'' '... % Excludes eADC
        'order by SeriesDate ASC']);  
series_eADC = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%Exponential Apparent Diffusion Coefficient%'' '... % Includes eADC
        'order by SeriesDate ASC']);   
series_FA = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%Fractional Aniso%'' '... % Includes FA
        'order by SeriesDate ASC']);   
series_AvgDC = mksqlite(['select * from Series where '...
        'SeriesDescription like ''%Average DC%'' '...  % Includes AvgDC
        'order by SeriesDate ASC']);   

%% DSC
series_rBV = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%rBV%'' '... % Includes rBV and rBV_corr
    'and SeriesDescription like ''%not%'' '... % Includes rBV
    'order by SeriesDate ASC']);
series_rBV_corr = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%rBV%'' '... % Includes rBV and rBV_corr
    'and SeriesDescription not like ''%not%'' '... % Excludes rBV
    'order by SeriesDate ASC']);
series_rBF = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%rBF%'' '... % Includes rBF 
    'order by SeriesDate ASC']);
series_MTT = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%MTT%'' '... % Includes MTT
    'order by SeriesDate ASC']);
series_Delay = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%Delay%'' '... % Includes Delay
    'order by SeriesDate ASC']);
series_K2 = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%Leakage (K2)%'' '... % Includes K2
    'order by SeriesDate ASC']);

%% DCE
series_Ktrans = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%K12%'' '... % Includes Ktrans
    'order by SeriesDate ASC,SeriesDescription DESC']);
series_Kep = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%K21%'' '... % Includes Kep
    'order by SeriesDate ASC,SeriesDescription DESC']);
series_Vp = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%Vp%'' '... % Includes Vp
    'order by SeriesDate ASC,SeriesDescription DESC']);
series_Ve = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%Ve (distribution volume)%'' '... % Includes Ve
    'order by SeriesDate ASC,SeriesDescription DESC']);
series_Wash_in = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%Wash-in%'' '... % Includes Wash-in
    'order by SeriesDate ASC,SeriesDescription DESC']);
series_Wash_out = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%Wash-out%'' '... % Includes Wash-out
    'order by SeriesDate ASC,SeriesDescription DESC']);
series_TTP = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%Time to peak%'' '... % Includes TTP
    'order by SeriesDate ASC,SeriesDescription DESC']);
series_AUC = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%Area under curve%'' '... % Includes AUC
    'order by SeriesDate ASC,SeriesDescription DESC']);
series_Peak = mksqlite(['select * from Series where '...
    'SeriesDescription like ''%Peak enhancement%'' '... % Includes Peak enhancement
    'order by SeriesDate ASC,SeriesDescription DESC']);

%% 
tic
x_step = .02;%abs(srow_x(srow_x(1:3)~=0)); % abs val of non-zero elements
y_step = .02;%abs(srow_y(srow_y(1:3)~=0));% mm
z_step = .02;%abs(srow_z(srow_z(1:3)~=0));
side = 1;%12; % side length in mm 
% Generate World Distances
[x_all y_all z_all] = meshgrid(...
    point(1)-(side/2):x_step:point(1)+(side/2),...
    point(2)-(side/2):y_step:point(2)+(side/2),...
    point(3)-(side/2):z_step:point(3)+(side/2));
x_diff = x_all - point(1);
y_diff = y_all - point(2);
z_diff = z_all - point(3);
xyz_dist = sqrt(x_diff.^2+y_diff.^2+z_diff.^2); % .02 Vox info  

load_point_info % loads point_info
for ii = 1%:size(point_info,1)
    point = cell2mat(point_info(ii,4:6));
    vcp = cell2mat(point_info(ii,7));

    if ~isequal('Forceps',vcp)% if not forceps
        % Random World -> .02 Vox
        mask_cyl = Mask_Cyl(ii, point_info, x_all);

        % .02 Vox -> .02 World
        idx_cyl = find(mask_cyl);
        x_cf = x_all(idx_cyl); 
        y_cf = y_all(idx_cyl);
        z_cf = z_all(idx_cyl);
    elseif isequal('Forceps',vcp) % 'Forceps'
        %% Forceps
        idx_forc = find(xyz_dist(:) < .5); % forceps .5 mm radius -> .02 Vox
        x_cf = x_all(idx_forc); % .02 Vox --> .02 World
        y_cf = y_all(idx_forc); 
        z_cf = z_all(idx_forc);
    end

    %% Sphere
    idx_sph = find(xyz_dist(:) < 2.5); % sphere 2.5 mm radius -> .02 Vox
    x_sph = x_all(idx_sph); % .02 Vox --> .02 World
    y_sph = y_all(idx_sph); 
    z_sph = z_all(idx_sph); 

    
    %% sort into this order
    % 5. T2FLAIR, - treat this separately
    prefix = script01_prefix;
    [cf_orig5 cf_high5 sph_orig5 sph_high5] = VOI_sample1(prefix,series_T2FLAIR,    ii,point_info, x_cf_FL,y_cf_FL,z_cf_FL, x_sph_FL,y_sph_FL,z_sph_FL);
    
    %% The rest
    prefix = script01_prefix;
    [cf_orig1 cf_high1 sph_orig1 sph_high1] = VOI_sample1(prefix,series_T1,         ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig2 cf_high2 sph_orig2 sph_high2] = VOI_sample1(prefix,series_T1post,     ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig3 cf_high3 sph_orig3 sph_high3] = VOI_sample1(prefix,series_T2,         ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig4 cf_high4 sph_orig4 sph_high4] = VOI_sample1(prefix,series_T2star,     ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig6 cf_high6 sph_orig6 sph_high6] = VOI_sample1(prefix,series_SWAN,       ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    
    prefix = script02_prefix;
    [cf_orig7 cf_high7 sph_orig7 sph_high7]  = VOI_sample2(prefix,series_ADC,   ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig8 cf_high8 sph_orig8 sph_high8]  = VOI_sample2(prefix,series_eADC,  ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig9 cf_high9 sph_orig9 sph_high9]  = VOI_sample2(prefix,series_FA,    ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig10 cf_high10 sph_orig10 sph_high10] = VOI_sample2(prefix,series_AvgDC, ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    
    prefix = script03_prefix;
    [cf_orig11 cf_high11 sph_orig11 sph_high11] = VOI_sample3(prefix,series_rBV,       ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig12 cf_high12 sph_orig12 sph_high12] = VOI_sample3(prefix,series_rBV_corr,  ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig13 cf_high13 sph_orig13 sph_high13] = VOI_sample3(prefix,series_rBF,       ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig14 cf_high14 sph_orig14 sph_high14] = VOI_sample3(prefix,series_MTT,       ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig15 cf_high15 sph_orig15 sph_high15] = VOI_sample3(prefix,series_Delay,     ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig16 cf_high16 sph_orig16 sph_high16] = VOI_sample3(prefix,series_K2,        ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    
    prefix = script03_prefix;
    [cf_orig17 cf_high17 sph_orig17 sph_high17] = VOI_sample3(prefix,series_Ktrans,           ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig18 cf_high18 sph_orig18 sph_high18] = VOI_sample3(prefix,series_Kep,              ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig19 cf_high19 sph_orig19 sph_high19] = VOI_sample3(prefix,series_Vp,               ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig20 cf_high20 sph_orig20 sph_high20] = VOI_sample3(prefix,series_Ve,               ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig21 cf_high21 sph_orig21 sph_high21] = VOI_sample3(prefix,series_Wash_in,          ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig22 cf_high22 sph_orig22 sph_high22] = VOI_sample3(prefix,series_Wash_out,         ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig23 cf_high23 sph_orig23 sph_high23] = VOI_sample3(prefix,series_TTP,              ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig24 cf_high24 sph_orig24 sph_high24] = VOI_sample3(prefix,series_AUC,              ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig25 cf_high25 sph_orig25 sph_high25] = VOI_sample3(prefix,series_Peak,             ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    

% Create cell arrays
voi_cf_orig = [...
    cf_orig1;   cf_orig2;   cf_orig3;   cf_orig4;   cf_orig5;   cf_orig6;...
    cf_orig7;   cf_orig8;   cf_orig9;   cf_orig10;  cf_orig11;  cf_orig12;...
    cf_orig13;  cf_orig14;  cf_orig15;  cf_orig16;  cf_orig17;  cf_orig18;...
    cf_orig19;  cf_orig20;  cf_orig21;  cf_orig22;  cf_orig23;  cf_orig24;  cf_orig25];

voi_cf_high = [...
    cf_high1;   cf_high2;   cf_high3;   cf_high4;   cf_high5;   cf_high6;...
    cf_high7;   cf_high8;   cf_high9;   cf_high10;  cf_high11;  cf_high12;...
    cf_high13;  cf_high14;  cf_high15;  cf_high16;  cf_high17;  cf_high18;...
    cf_high19;  cf_high20;  cf_high21;  cf_high22;  cf_high23;  cf_high24;  cf_high25];

voi_sph_orig = [...
    sph_orig1;   sph_orig2;   sph_orig3;   sph_orig4;   sph_orig5;   sph_orig6;...
    sph_orig7;   sph_orig8;   sph_orig9;   sph_orig10;  sph_orig11;  sph_orig12;...
    sph_orig13;  sph_orig14;  sph_orig15;  sph_orig16;  sph_orig17;  sph_orig18;...
    sph_orig19;  sph_orig20;  sph_orig21;  sph_orig22;  sph_orig23;  sph_orig24;  sph_orig25];

voi_sph_high = [...
    sph_high1;   sph_high2;   sph_high3;   sph_high4;   sph_high5;   sph_high6;...
    sph_high7;   sph_high8;   sph_high9;   sph_high10;  sph_high11;  sph_high12;...
    sph_high13;  sph_high14;  sph_high15;  sph_high16;  sph_high17;  sph_high18;...
    sph_high19;  sph_high20;  sph_high21;  sph_high22;  sph_high23;  sph_high24;  sph_high25];


    %% V. Use VOI mask to save mean, sd, histogram in .mat file
    ptno = point_info(ii,1);
    ptno_site = point_info(ii,2);
    script04_results = [script04_prefix 'Results/'];
    if ~isdir(script04_results)
        mkdir(script04_results)
    end

    % save .mat files in 1 directory
    if ~isequal('Forceps',vcp)% if not forceps
        save([script04_results ptno_site '_CO.mat'], 'voi_cf_orig') 
        save([script04_results ptno_site '_CH.mat'], 'voi_cf_high') 
    elseif isequal('Forceps',vcp) % 'Forceps'
        save([script04_results ptno_site '_FO.mat'], 'voi_cf_orig') 
        save([script04_results ptno_site '_FH.mat'], 'voi_cf_high') 
    end
    save([script04_results ptno_site '_SO.mat'], 'voi_sph_orig') 
    save([script04_results ptno_site '_SH.mat'], 'voi_sph_high') 

end

time_hrs=toc/60/60 % 26.74 min




