%% I. Load stuff
clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql

%% 
tic
x_step = .02; 
y_step = x_step;
z_step = x_step;
side = 12; % side length in mm 
% Generate World Distances

load_point_info % loads point_info
for ii = 1:size(point_info,1)
    [ii point_info(ii,2)]
    point = cell2mat(point_info(ii,4:6));
    [x_all y_all z_all] = meshgrid(...
        point(1)-(side/2):x_step:point(1)+(side/2),...
        point(2)-(side/2):y_step:point(2)+(side/2),...
        point(3)-(side/2):z_step:point(3)+(side/2));
    x_diff = x_all - point(1);
    y_diff = y_all - point(2);
    z_diff = z_all - point(3);
    xyz_dist = sqrt(x_diff.^2+y_diff.^2+z_diff.^2); % .02 Vox info  

    vcp = cell2mat(point_info(ii,7));
    if ~isequal('Forceps',vcp)% if not forceps
        % Random World -> .02 Vox
        mask_cyl = Mask_Cyl(ii, point_info, series_T2_register, x_all,x_step,y_step,z_step,side);

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

    display('Coordinates done')
    
    %% sort into this order
    % 5. T2FLAIR, - treat this separately
    %tic
    lot=1;
    prefix = script01_prefix;
    x_cf_FL = x_cf;
    y_cf_FL = y_cf;
    z_cf_FL = z_cf;
    x_sph_FL = x_sph;
    y_sph_FL = y_sph;
    z_sph_FL = z_sph;
    [cf_orig5 cf_high5 sph_orig5 sph_high5] = VOI_sample(1  ,prefix,series_T2FLAIR,     ii,point_info, x_cf_FL,y_cf_FL,z_cf_FL, x_sph_FL,y_sph_FL,z_sph_FL);
    [cf_orig3 cf_high3 sph_orig3 sph_high3] = VOI_sample(0  ,prefix,series_T2_register,ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    
    %% The rest
    lot=1;
    prefix = script01_prefix;
    [cf_orig1 cf_high1 sph_orig1 sph_high1] = VOI_sample(lot,prefix,series_T1,         ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig2 cf_high2 sph_orig2 sph_high2] = VOI_sample(lot,prefix,series_T1post,     ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig4 cf_high4 sph_orig4 sph_high4] = VOI_sample(lot,prefix,series_T2star,     ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig6 cf_high6 sph_orig6 sph_high6] = VOI_sample(lot,prefix,series_SWAN,       ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    display('T1 T2 Swan done')
    
    lot=2;
    prefix = script02_prefix;
    [cf_orig7 cf_high7 sph_orig7 sph_high7]  = VOI_sample(lot,prefix,series_ADC,   ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig8 cf_high8 sph_orig8 sph_high8]  = VOI_sample(lot,prefix,series_eADC,  ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig9 cf_high9 sph_orig9 sph_high9]  = VOI_sample(lot,prefix,series_FA,    ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig10 cf_high10 sph_orig10 sph_high10] = VOI_sample(lot,prefix,series_AvgDC, ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    display('DWI DTI done')
    
    lot=3;
    prefix = script03_prefix;
    [cf_orig11 cf_high11 sph_orig11 sph_high11] = VOI_sample(lot,prefix,series_rBV,       ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig12 cf_high12 sph_orig12 sph_high12] = VOI_sample(lot,prefix,series_rBV_corr,  ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig13 cf_high13 sph_orig13 sph_high13] = VOI_sample(lot,prefix,series_rBF,       ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig14 cf_high14 sph_orig14 sph_high14] = VOI_sample(lot,prefix,series_MTT,       ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig15 cf_high15 sph_orig15 sph_high15] = VOI_sample(lot,prefix,series_Delay,     ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig16 cf_high16 sph_orig16 sph_high16] = VOI_sample(lot,prefix,series_K2,        ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    display('DSC Done')
    
    lot=3;
    prefix = script03_prefix;
    [cf_orig17 cf_high17 sph_orig17 sph_high17] = VOI_sample(lot,prefix,series_Ktrans,           ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig18 cf_high18 sph_orig18 sph_high18] = VOI_sample(lot,prefix,series_Kep,              ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig19 cf_high19 sph_orig19 sph_high19] = VOI_sample(lot,prefix,series_Vp,               ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig20 cf_high20 sph_orig20 sph_high20] = VOI_sample(lot,prefix,series_Ve,               ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig21 cf_high21 sph_orig21 sph_high21] = VOI_sample(lot,prefix,series_Wash_in,          ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig22 cf_high22 sph_orig22 sph_high22] = VOI_sample(lot,prefix,series_Wash_out,         ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig23 cf_high23 sph_orig23 sph_high23] = VOI_sample(lot,prefix,series_TTP,              ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig24 cf_high24 sph_orig24 sph_high24] = VOI_sample(lot,prefix,series_AUC,              ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    [cf_orig25 cf_high25 sph_orig25 sph_high25] = VOI_sample(lot,prefix,series_Peak,             ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph);
    display('DCE done')
    %time_min = toc/60 % 5.26 min!


% Create cell arrays
voi_cf_orig = [...
    cf_orig1;   cf_orig2;   cf_orig3;   cf_orig4;   cf_orig5;   cf_orig6;... %t1, t2, swan
    cf_orig7;   cf_orig8;   cf_orig9;   cf_orig10;... % adc, fa
    cf_orig11;  cf_orig12;  cf_orig13;  cf_orig14;  cf_orig15;  cf_orig16;... %dsc
    cf_orig17;  cf_orig18;  cf_orig19;  cf_orig20;  cf_orig21;  cf_orig22;... %dce
    cf_orig23;  cf_orig24;  cf_orig25];

voi_cf_high = [...
    cf_high1;   cf_high2;   cf_high3;   cf_high4;   cf_high5;   cf_high6;...
    cf_high7;   cf_high8;   cf_high9;   cf_high10;...
    cf_high11;  cf_high12;  cf_high13;  cf_high14;  cf_high15;  cf_high16;...
    cf_high17;  cf_high18;  cf_high19;  cf_high20;  cf_high21;  cf_high22;...
    cf_high23;  cf_high24;  cf_high25];

voi_sph_orig = [...
    sph_orig1;   sph_orig2;   sph_orig3;   sph_orig4;   sph_orig5;   sph_orig6;...
    sph_orig7;   sph_orig8;   sph_orig9;   sph_orig10;...
    sph_orig11;  sph_orig12;  sph_orig13;  sph_orig14;  sph_orig15;  sph_orig16;...
    sph_orig17;  sph_orig18;  sph_orig19;  sph_orig20;  sph_orig21;  sph_orig22;...
    sph_orig23;  sph_orig24;  sph_orig25];

voi_sph_high = [...
    sph_high1;   sph_high2;   sph_high3;   sph_high4;   sph_high5;   sph_high6;...
    sph_high7;   sph_high8;   sph_high9;   sph_high10;...
    sph_high11;  sph_high12;  sph_high13;  sph_high14;  sph_high15;  sph_high16;...
    sph_high17;  sph_high18;  sph_high19;  sph_high20;  sph_high21;  sph_high22;...
    sph_high23;  sph_high24;  sph_high25];


    %% V. Use VOI mask to save mean, sd, histogram in .mat file
    ptno = cell2mat(point_info(ii,1));
    ptno_site = cell2mat(point_info(ii,2));
    
    script04_results = [script04_prefix 'Results/'];
    script04_results_ptno = [script04_results ptno '/'];
    script04_results_ptno_site = [script04_results_ptno ptno_site '/'];
   
    if ~isdir(script04_results)
        mkdir(script04_results)
    end
    if ~isdir(script04_results_ptno)
        mkdir(script04_results_ptno)
    end 
    if ~isdir(script04_results_ptno_site)
        mkdir(script04_results_ptno_site)
    end 
    
    %% Create mask to Check coordinates - make sure they're in the right place
    jj=str2num(cell2mat(point_info(ii,1)));
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_T2_register(jj));

    % Original Rez
    file = load_nii([ script01_prefix '01_N4/N4_' ptno_Descrip_UID '.nii.gz' ]);
    srow_x = file.hdr.hist.srow_x; % [ITK k]
    srow_y = file.hdr.hist.srow_y; % [ITK i]
    srow_z = file.hdr.hist.srow_z; % [ITK j]
    mask_cf = zeros(size(file.img)); % Create masks
    mask_sph = zeros(size(file.img));
   
    % For original rez - Find coordinates, convert to single index, put 1
    [i_cf j_cf k_cf] = World_to_Vox(x_cf, y_cf, z_cf, srow_x, srow_y, srow_z); % .02 World -> .5469 Vox 
    [i_sph j_sph k_sph] = World_to_Vox(x_sph, y_sph, z_sph, srow_x, srow_y, srow_z);% .02 World --> .5469 Vox
    cf_sample = sub2ind(size(file.img),i_cf,j_cf,k_cf); % .5469 Vox -> Values
    sph_sample = sub2ind(size(file.img),i_sph,j_sph,k_sph); % .5469 Vox -> Values
    mask_cf(cf_sample) = 1;
    mask_sph(sph_sample) = 1; 

    % Display Rez
    Generate_Display
    file_d = load_nii([ script01_prefix '11_Resample/Resample_' ptno_site '_' ptno_Descrip_UID '.nii.gz' ]);
    srow_x_d = file_d.hdr.hist.srow_x; % [ITK k]
    srow_y_d = file_d.hdr.hist.srow_y; % [ITK i]
    srow_z_d = file_d.hdr.hist.srow_z; % [ITK j]
    mask_cf_d = zeros(size(file_d.img));
    mask_sph_d = zeros(size(file_d.img));
   
    % For display rez - Find coordinates, convert to single index, put 1
    [i_cf_d j_cf_d k_cf_d] = World_to_Vox(x_cf, y_cf, z_cf, srow_x_d, srow_y_d, srow_z_d); % .02 World -> .5469 Vox 
    [i_sph_d j_sph_d k_sph_d] = World_to_Vox(x_sph, y_sph, z_sph, srow_x_d, srow_y_d, srow_z_d);% .02 World --> .5469 Vox
    cf_sample_d = sub2ind(size(file_d.img),i_cf_d,j_cf_d,k_cf_d); % .5469 Vox -> Values
    sph_sample_d = sub2ind(size(file_d.img),i_sph_d,j_sph_d,k_sph_d); % .5469 Vox -> Values
    mask_cf_d(cf_sample_d) = 1;
    mask_sph_d(sph_sample_d) = 1; 
    display('Masks and Display volume done')

    % save .mat files in 1 directory
    if ~isequal('Forceps',vcp)% if not forceps
        
        save([script04_results ptno_site '_CO.mat'], 'voi_cf_orig') 
        save([script04_results ptno_site '_CH.mat'], 'voi_cf_high') 
        
        holder_cf = make_nii(mask_cf,size(file.img),[srow_x(4) srow_y(4) srow_z(4)],64,['mask_' ptno_site '_CO']);
        holder_cf.hdr.dime.pixdim = file.hdr.dime.pixdim;
        holder_cf.hdr.hist = file.hdr.hist;
        save_nii(holder_cf,[script04_results 'Mask_' ptno_site '_CO.nii.gz']) 
        
        holder_cf_d = make_nii(mask_cf_d,size(file_d.img),[srow_x_d(4) srow_y_d(4) srow_z_d(4)],64,['mask_' ptno_site '_CH']);
        holder_cf_d.hdr.dime.pixdim = file_d.hdr.dime.pixdim;
        holder_cf_d.hdr.hist = file_d.hdr.hist;
        save_nii(holder_cf_d,[script04_results 'Mask_' ptno_site '_CH.nii.gz']) 
        
    elseif isequal('Forceps',vcp) % 'Forceps'
        
        save([script04_results ptno_site '_FO.mat'], 'voi_cf_orig') 
        save([script04_results ptno_site '_FH.mat'], 'voi_cf_high') 
        
        holder_cf = make_nii(mask_cf,size(file.img),[srow_x(4) srow_y(4) srow_z(4)],64,['mask_' ptno_site '_FO']);
        holder_cf.hdr.dime.pixdim = file.hdr.dime.pixdim;
        holder_cf.hdr.hist = file.hdr.hist;
        save_nii(holder_cf,[script04_results 'Mask_' ptno_site '_FO.nii.gz']) 
        
        holder_cf_d = make_nii(mask_cf_d,size(file_d.img),[srow_x_d(4) srow_y_d(4) srow_z_d(4)],64,['mask_' ptno_site '_FH']);
        holder_cf_d.hdr.dime.pixdim = file_d.hdr.dime.pixdim;
        holder_cf_d.hdr.hist = file_d.hdr.hist;
        save_nii(holder_cf_d,[script04_results 'Mask_' ptno_site '_FH.nii.gz']) 
        
    end
    
    save([script04_results ptno_site '_SO.mat'], 'voi_sph_orig') 
    save([script04_results ptno_site '_SH.mat'], 'voi_sph_high') 
    
    holder_sph = make_nii(mask_sph,size(file.img),[srow_x(4) srow_y(4) srow_z(4)],64,['mask_' ptno_site '_SO']);
    holder_sph.hdr.dime.pixdim = file.hdr.dime.pixdim;
    holder_sph.hdr.hist = file.hdr.hist;
    save_nii(holder_sph,[script04_results 'Mask_' ptno_site '_SO.nii.gz']) 

    holder_sph_d = make_nii(mask_sph_d,size(file_d.img),[srow_x_d(4) srow_y_d(4) srow_z_d(4)],64,['mask_' ptno_site '_SH']);
    holder_sph_d.hdr.dime.pixdim = file_d.hdr.dime.pixdim;
    holder_sph_d.hdr.hist = file_d.hdr.hist;
    save_nii(holder_sph_d,[script04_results 'Mask_' ptno_site '_SH.nii.gz']) 
    
    
end

toc/60 % 114.47 min for 14 points
toc/60/60 % 1.91 hrs for 14 points




