function [cf_orig cf_high sph_orig sph_high]= VOI_sample(lot, prefix, series_sql,ii,point_info, x_cf,y_cf,z_cf, x_sph,y_sph,z_sph)

path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql

if lot==0 %t2 register - will I ever use this one?
    img_dir = [ prefix '01_N4/N4_'];
    % img_dir = [ prefix '09_GaussNorm/GN_' ];
elseif lot==1 % t2 histmatch, t1, swan
    img_dir = [ prefix '07_Warped/Warped_'  ];
    % img_dir = [ prefix '09_GaussNorm/GN_' ];
elseif lot==2 % adc,fa
    img_dir = [ prefix '07_Warped/Warped_'  ];
    % img_dir = [ prefix '09_GaussNorm/GN_' ];
elseif lot==3 % dce/dsc
    img_dir = [prefix '08_Warped/Warped_'];
    % img_dir = [ prefix '09_GaussNorm/GN_' ];
end

ptno = cell2mat(point_info(ii,1));
ptno_site = cell2mat(point_info(ii,2));
aim = cell2mat(point_info(ii,3));
vcp = cell2mat(point_info(ii,7)); % vector control point
    
script04_results_ptno = [script04_prefix 'Results/' ptno '/'];
script04_results_ptno_site = [script04_results_ptno ptno_site '/'];
if ~isdir(script04_results_ptno)
    mkdir(script04_results_ptno)
end 
if ~isdir(script04_results_ptno_site)
    mkdir(script04_results_ptno_site)
end 

%% ptno --> Correct sql entry
% Ptno -> StudyInstanceUID -> series_sql entry
% series_sql = series_T2_Register; % for troubleshooting
[field study]= find(cell2mat(cellfun(@(x) isequal(studies_all(str2num(ptno)).StudyInstanceUID,x),struct2cell(series_sql),'UniformOutput',false)));

for jj = study'
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(jj));
    
    if lot==0 || lot==1
        conv_factor = double(1);
    elseif lot==2 % adc,fa
        conv_factor = double(1e-6);
    elseif lot==3 % dce/dsc
        images = mksqlite(['select * from Images where SeriesInstanceUID = ''' UID ''' order by Filename ASC' ]);
        tmp2=dicominfo(images(1).Filename,'dictionary','gems-dicom-dict.txt');
        conv_factor = double(tmp2.Private_0077_1001);
    end
    
    file = load_nii([ img_dir ptno_Descrip_UID '.nii.gz' ]);
    srow_x = file.hdr.hist.srow_x; % [ITK k]
    srow_y = file.hdr.hist.srow_y; % [ITK i]
    srow_z = file.hdr.hist.srow_z; % [ITK j]

    % Cylinder/Forceps
    [i_cf, j_cf, k_cf] = World_to_Vox(x_cf, y_cf, z_cf, srow_x, srow_y, srow_z); % .02 World -> .5469 Vox 
    cf_sample = sub2ind(size(file.img),i_cf,j_cf,k_cf); % .5469 Vox -> Values
    cf_orig_values = double(file.img(unique(cf_sample)))*conv_factor; % Original Rez
    cf_high_values = double(file.img(cf_sample))*conv_factor; % High Rez
    [cf_orig_mean, cf_orig_std, cf_orig_min, cf_orig_max, cf_orig_range, cf_orig_numel] = Generate_Stats(cf_orig_values);
    [cf_high_mean, cf_high_std, cf_high_min, cf_high_max, cf_high_range, cf_high_numel] = Generate_Stats(cf_high_values);

    % Sphere
    [i_sph, j_sph, k_sph] = World_to_Vox(x_sph, y_sph, z_sph, srow_x, srow_y, srow_z);% .02 World --> .5469 Vox
    sph_sample = sub2ind(size(file.img),i_sph,j_sph,k_sph); % .5469 Vox -> Values
    sph_orig_values = double(file.img(unique(sph_sample)))*conv_factor; % Original Rez
    sph_high_values = double(file.img(sph_sample))*conv_factor; % High Rez
    [sph_orig_mean, sph_orig_std, sph_orig_min, sph_orig_max, sph_orig_range, sph_orig_numel] = Generate_Stats(sph_orig_values);
    [sph_high_mean, sph_high_std, sph_high_min, sph_high_max, sph_high_range, sph_high_numel] = Generate_Stats(sph_high_values);


    %% Generate Output values
    if ~isequal(vcp,'Forceps')
        ptno_site_shape1 = [ptno_site '_CO'];
        ptno_site_shape2 = [ptno_site '_CH'];
    elseif isequal(vcp,'Forceps')
        ptno_site_shape1 = [ptno_site '_FO'];
        ptno_site_shape2 = [ptno_site '_FH'];
    end
    ptno_site_shape3 = [ptno_site '_SO'];
    ptno_site_shape4 = [ptno_site '_SH'];
    
    if ~exist('cf_orig')
        cf_orig = {ptno ptno_site_shape1 aim Descrip...
            cf_orig_mean cf_orig_std...
            cf_orig_min cf_orig_max...
            cf_orig_range cf_orig_numel};
        cf_high = {ptno ptno_site_shape2 aim Descrip...
            cf_high_mean cf_high_std...
            cf_high_min cf_high_max...
            cf_high_range cf_high_numel}; 
        sph_orig = {ptno ptno_site_shape3 aim Descrip...
            sph_orig_mean sph_orig_std...
            sph_orig_min sph_orig_max...
            sph_orig_range sph_orig_numel};
        sph_high = {ptno ptno_site_shape4 aim Descrip...
            sph_high_mean sph_high_std...
            sph_high_min sph_high_max...
            sph_high_range sph_high_numel}; 
    else % if it's already there - Ktrans etc where there's 4 versions of the same map
        cf_orig = [cf_orig; {ptno ptno_site_shape1 aim Descrip...
            cf_orig_mean cf_orig_std...
            cf_orig_min cf_orig_max...
            cf_orig_range cf_orig_numel}];
        cf_high = [cf_high; {ptno ptno_site_shape2 aim Descrip...
            cf_high_mean cf_high_std...
            cf_high_min cf_high_max...
            cf_high_range cf_high_numel}];
        sph_orig = [sph_orig; {ptno ptno_site_shape3 aim Descrip...
            sph_orig_mean sph_orig_std...
            sph_orig_min sph_orig_max...
            sph_orig_range sph_orig_numel}];
        sph_high = [sph_high; {ptno ptno_site_shape4 aim Descrip...
            sph_high_mean sph_high_std...
            sph_high_min sph_high_max...
            sph_high_range sph_high_numel}]; 
    end

    % Generate/Save Histograms
    label_x = 'Value';
    label_y = '# Voxels';
    if ~isequal('Forceps',vcp)% if not forceps
        figure;		hist(cf_orig_values);        xlabel(label_x);        ylabel(label_y)
        title({     [ ptno_site ', ' Descrip];...
                    ['mean = '      num2str(cf_orig_mean)...
                    ' +/- '         num2str(cf_orig_std)...
                    ', # vox = '    num2str(cf_orig_numel) ', CO']}, 'Interpreter' ,'none')
        saveas(gcf, [script04_results_ptno_site		ptno_site '_Histogram_' Descrip '_CO' ],'png')
        close all

        figure;		hist(cf_high_values);        xlabel(label_x);        ylabel(label_y)
        title({     [ ptno_site ', ' Descrip];...
                    ['mean = '      num2str(cf_high_mean)...
                    ' +/- '         num2str(cf_high_std)...
                    ', # vox = '    num2str(cf_high_numel) ', CH']},'Interpreter','none')
        saveas(gcf, [script04_results_ptno_site		ptno_site '_Histogram_' Descrip '_CH' ],'png')
        close all

    elseif isequal('Forceps',vcp) % 'Forceps'
        figure;		hist(cf_orig_values);        xlabel(label_x);        ylabel(label_y)
        title({     [ ptno_site ', ' Descrip];...
                    ['mean = '      num2str(cf_orig_mean)...
                    ' +/- '         num2str(cf_orig_std)...
                    ', # vox = '    num2str(cf_orig_numel) ', FO']},'Interpreter','none')
        saveas(gcf, [script04_results_ptno_site		ptno_site '_Histogram_' Descrip '_FO' ],'png')
        close all

        figure; 	hist(cf_high_values);        xlabel(label_x);        ylabel(label_y)
        title({     [ ptno_site ', ' Descrip];...
                    ['mean = '      num2str(cf_high_mean)...
                    ' +/- '         num2str(cf_high_std)...
                    ', # vox = '    num2str(cf_high_numel) ', FH']},'Interpreter','none')
        saveas(gcf, [script04_results_ptno_site		ptno_site '_Histogram_' Descrip '_FH' ],'png')
        close all
    end

    figure;	hist(sph_orig_values);    xlabel(label_x);    ylabel(label_y)
    title({     [ ptno_site ', ' Descrip];...
                ['mean = '      num2str(sph_orig_mean)...
                ' +/- '         num2str(sph_orig_std)...
                ', # vox = '    num2str(sph_orig_numel) ', SO']},'Interpreter','none')
    saveas(gcf,  [script04_results_ptno_site	ptno_site '_Histogram_' Descrip '_SO' ],'png')
    close all

    figure;	hist(sph_high_values);    xlabel(label_x);    ylabel(label_y)
    title({     [ ptno_site ', ' Descrip];...
                ['mean = '      num2str(sph_high_mean)...
                ' +/- '         num2str(sph_high_std)...
                ', # vox = '    num2str(sph_high_numel) ', SH']},'Interpreter','none')
    saveas(gcf,  [script04_results_ptno_site	ptno_site '_Histogram_' Descrip '_SH' ],'png')
    close all
end

end


