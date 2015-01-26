function [csf gm wm fat] = ROI_sample(lot, prefix, series_sql, ROI_ptno, idx1, idx2, idx3, idx4) 

path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql
script04_results_ptno = [script04_prefix 'Results/' ROI_ptno '/'];
if ~isdir(script04_results_ptno)
    mkdir(script04_results_ptno)
end 

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

[field, study]= find(cell2mat(cellfun(@(x) isequal(studies_all(str2num(ROI_ptno)).StudyInstanceUID,x),struct2cell(series_sql),'UniformOutput',false)));
for jj = study'
    [img_ptno, img_Descrip, img_UID, img_ptno_Descrip_UID] = Generate_Label(series_sql(jj));

    if lot==0 || lot==1
        conv_factor = double(1);
    elseif lot==2 % adc,fa
        conv_factor = double(1e-6);
    elseif lot==3 % dce/dsc
        images = mksqlite(['select * from Images where SeriesInstanceUID = ''' img_UID ''' order by Filename ASC' ]);
        tmp2=dicominfo(images(1).Filename,'dictionary','gems-dicom-dict.txt');
        conv_factor = double(tmp2.Private_0077_1001);
    end

    file = load_nii([ img_dir  img_ptno_Descrip_UID '.nii.gz' ]);
    csf_values = double(file.img(idx1))*conv_factor;
    gm_values = double(file.img(idx2))*conv_factor;
    wm_values = double(file.img(idx3))*conv_factor;
    fat_values = double(file.img(idx4))*conv_factor;
    [csf_mean csf_std csf_min csf_max csf_range csf_numel] = Generate_Stats(csf_values);
    [gm_mean  gm_std  gm_min  gm_max  gm_range  gm_numel] = Generate_Stats(gm_values);
    [wm_mean  wm_std  wm_min  wm_max  wm_range  wm_numel] = Generate_Stats(wm_values);
    [fat_mean fat_std fat_min fat_max fat_range fat_numel] = Generate_Stats(fat_values);

    if ~exist('csf')
        csf = {img_ptno img_Descrip 'CSF' csf_mean csf_std csf_min csf_max csf_range csf_numel};
        gm = {img_ptno img_Descrip 'GM' gm_mean  gm_std  gm_min  gm_max  gm_range  gm_numel};
        wm = {img_ptno img_Descrip 'WM' wm_mean  wm_std  wm_min  wm_max  wm_range  wm_numel};
        fat = {img_ptno img_Descrip 'Fat' fat_mean fat_std fat_min fat_max fat_range fat_numel};
    else
        csf = [csf;{ptno Descrip 'CSF' csf_mean csf_std csf_min csf_max csf_range csf_numel}];
        gm = [gm;{ptno Descrip 'GM' gm_mean  gm_std  gm_min  gm_max  gm_range  gm_numel}];
        wm = [wm;{ptno Descrip 'WM' wm_mean  wm_std  wm_min  wm_max  wm_range  wm_numel}];
        fat = [fat;{ptno Descrip 'Fat' fat_mean fat_std fat_min fat_max fat_range fat_numel}];
    end

    %% Generate/Save Histogram - get this done after grant
    label_x = 'Value';
    label_y = '# Voxels';

    figure;     hist(csf_values);       xlabel(label_x);        ylabel(label_y)
    title({     [ img_ptno ', ' img_Descrip ', CSF' ];...
                ['mean = '      num2str(csf_mean)...
                ' +/- '         num2str(csf_std)...
                ', # vox = '    num2str(csf_numel)]}, 'Interpreter' ,'none')
    saveas(gcf, [script04_results_ptno      'Histogram_' img_ptno '_' img_Descrip '_CSF' ],'png')
    close all

    figure;     hist(gm_values);        xlabel(label_x);        ylabel(label_y)
    title({     [ img_ptno ', ' img_Descrip ', GM' ];...
                ['mean = '      num2str(gm_mean)...
                ' +/- '         num2str(gm_std)...
                ', # vox = '    num2str(gm_numel)]},'Interpreter','none')
    saveas(gcf, [script04_results_ptno      'Histogram_' img_ptno '_' img_Descrip '_GM' ],'png')
    close all

    figure;     hist(wm_values);        xlabel(label_x);        ylabel(label_y)
    title({     [ img_ptno ', ' img_Descrip ', WM' ];...
                ['mean = '      num2str(wm_mean)...
                ' +/- '         num2str(wm_std)...
                ', # vox = '    num2str(wm_numel)]},'Interpreter','none')
    saveas(gcf,  [script04_results_ptno     'Histogram_' img_ptno '_' img_Descrip '_WM' ],'png')
    close all

    figure;     hist(fat_values);       xlabel(label_x);        ylabel(label_y)
    title({     [ img_ptno ', ' img_Descrip ', Fat' ];...
                ['mean = '      num2str(fat_mean)...
                ' +/- '         num2str(fat_std)...
                ', # vox = '    num2str(fat_numel)]},'Interpreter','none')
    saveas(gcf,  [script04_results_ptno     'Histogram_' img_ptno '_' img_Descrip '_Fat' ],'png')
    close all
end


end
                
                
                