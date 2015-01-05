function [cf_orig cf_high sph_orig sph_high]= VOI_sample1(prefix, series_sql,ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph)

mksqlite('open', 'ctkDICOM.sql' );

% Cell output possible?
%% Load
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';
ptno = point_info(ii,1);
ptno_site = point_info(ii,2);
script04_results_ptno_site = [script04_prefix 'Results/' ptno_site '/'];
if ~isdir(script04_results_ptno_site)
    mkdir(script04_results_ptno_site)
end 

%% ptno --> Correct sql entry
% ptno -> MRN 
mrn = mrn_to_ptno(ptno);
% MRN/PatientID -> series_sql entry
[field study]= find(cell2mat(cellfun(@(x) isequal(mrn,x),struct2cell(series_sql),'UniformOutput',false)));

for jj = study'
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(jj));
    images = mksqlite(['select * from Images where SeriesInstanceUID = ''' UID ''' order by Filename ASC' ]);
    tmp2=dicominfo(images(1).Filename,'dictionary','gems-dicom-dict.txt');

    file = load_nii([ prefix '04_Warped/Warped_' ptno_Descrip_UID '.nii.gz' ]);
    % file = load_nii(['test.nii.gz']);
    srow_x = file.hdr.hist.srow_x; % [ITK k]
    srow_y = file.hdr.hist.srow_y; % [ITK i]
    srow_z = file.hdr.hist.srow_z; % [ITK j]


    % Cylinder/Forceps
    [i_cyl j_cyl k_cyl] = World_to_Vox(x_cf, y_cf, z_cf, srow_x, srow_y, srow_z); % .02 World -> .5469 Vox 
    cf_sample = sub2ind(size(file.img),i_cyl,j_cyl,k_cyl); % .5469 Vox -> Values
    cf_orig = file.img(unique(cf_sample))*tmp2.Private_0077_1001; % Original Rez
    cf_orig_mean = mean(cf_orig,1);
    cf_orig_std = std(cf_orig,0,1);
    cf_orig_min = min(cf_orig,[],1);
    cf_orig_max = max(cf_orig,[],1);
    cf_orig_range = cf_orig_max - cf_orig_min;
    cf_orig_numel = size(cf_orig,1);  
    cf_high = file.img(cf_sample)*tmp2.Private_0077_1001; % High Rez
    cf_high_mean = mean(cf_high,1);
    cf_high_std = std(cf_high,0,1);
    cf_high_min = min(cf_high,[],1);
    cf_high_max = max(cf_high,[],1);
    cf_high_range = cf_high_max - cf_high_min;
    cf_high_numel = size(cf_high,1);

    % Sphere
    [i_sph j_sph k_sph] = World_to_Vox(x_sph, y_sph, z_sph, srow_x, srow_y, srow_z);% .02 World --> .5469 Vox
    sph_sample = sub2ind(size(file.img),i_sph,j_sph,k_sph); % .5469 Vox -> Values
    sph_orig = file.img(unique(sph_sample))*tmp2.Private_0077_1001; % Original Rez
    sph_orig_mean = mean(sph_orig,1);
    sph_orig_std = std(sph_orig,0,1);
    sph_orig_min = min(sph_orig,[],1);
    sph_orig_max = max(sph_orig,[],1);
    sph_orig_range = sph_orig_max - sph_orig_min;
    sph_orig_numel = size(sph_orig,1);
    sph_high = file.img(sph_sample)*tmp2.Private_0077_1001; % High Rez
    sph_high_mean = mean(sph_high,1);
    sph_high_std = std(sph_high,0,1);
    sph_high_min = min(sph_high,[],1);
    sph_high_max = max(sph_high,[],1);
    sph_high_range = sph_high_max - sph_high_min;
    sph_high_numel = size(sph_high,1);


    %% Generate Output values
    cf_orig = {ptno ptno_site Descrip...
        cf_orig_mean cf_orig_std...
        cf_orig_min cf_orig_max...
        cf_orig_range cf_orig_numel};
    cf_high = {ptno ptno_site Descrip...
        cf_high_mean cf_high_std...
        cf_high_min cf_high_max...
        cf_high_range cf_high_numel};
    sph_orig = {ptno ptno_site Descrip...
        sph_orig_mean sph_orig_std...
        sph_orig_min sph_orig_max...
        sph_orig_range sph_orig_numel};
    sph_high = {ptno ptno_site Descrip...
        sph_high_mean sph_high_std...
        sph_high_min sph_high_max...
        sph_high_range sph_high_numel};

    %% Generate/Save Histograms
    vcp = cell2mat(point_info(ii,7));
    if ~isequal('Forceps',vcp)% if not forceps
        h1 = figure;
        hist(cf_orig)
        title([ ptno ' ' site ', ' Descrip...
            ', mean = ' num2str(cf_orig_mean)...
            ' +/- ' num2str(cf_orig_std)...
            ', # pix = ' num2str(cf_orig_numel) ' CO'])
        saveas(h1, [script04_results_ptno_site...
            ptno '_' site '_Histogram_' Descrip '_CO' ],'png')

        h2 = figure;
        hist(cf_high)
        title([ ptno ' ' site ', ' Descrip...
            ', mean = ' num2str(cf_high_mean)...
            ' +/- ' num2str(cf_high_std)...
            ', # pix = ' num2str(cf_high_numel) ' CH'])
        saveas(h2, [script04_results_ptno_site...
            ptno '_' site '_Histogram_' Descrip '_CH' ],'png')

    elseif isequal('Forceps',vcp) % 'Forceps'
        h1 = figure;
        hist(cf_orig)
        title([ ptno ' ' site ', ' Descrip...
            ', mean = ' num2str(cf_orig_mean)...
            ' +/- ' num2str(cf_orig_std)...
            ', # pix = ' num2str(cf_orig_numel) ' FO'])
        saveas(h1, [script04_results_ptno_site...
            ptno '_' site '_Histogram_' Descrip '_FO' ],'png')

        h2 = figure;
        hist(cf_high)
        title([ ptno ' ' site ', ' Descrip...
            ', mean = ' num2str(cf_high_mean)...
            ' +/- ' num2str(cf_high_std)...
            ', # pix = ' num2str(cf_high_numel) ' FH'])
        saveas(h2, [script04_results_ptno_site...
            ptno '_' site '_Histogram_' Descrip '_FH' ],'png')

    end

    h3 = figure;
    hist(sph_orig)
    title([ ptno_site ', ' Descrip...
        ', mean = ' num2str(sph_orig_mean)...
        ' +/- ' num2str(sph_orig_std)...
        ', # pix = ' num2str(sph_orig_numel) ' SO'])
    saveas(h3, [script04_results_ptno_site...
        ptno_site '_Histogram_' Descrip '_SO' ],'png')

    h4 = figure;
    hist(sph_high)
    title([ ptno_site ', ' Descrip...
        ', mean = ' num2str(sph_high_mean)...
        ' +/- ' num2str(sph_high_std)...
        ', # pix = ' num2str(sph_high_numel) ' SH'])
    saveas(h4, [script04_results_ptno_site...
        ptno_site '_Histogram_' Descrip '_SH' ],'png')

    close all
end

end
