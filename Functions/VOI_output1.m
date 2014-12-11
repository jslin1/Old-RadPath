
function  VOI_output1(temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)

    path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
    path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
    path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
    script01_prefix = 'Script01_T1_T2_SWAN/';
    script02_prefix = 'Script02_DWI_DTI/';
    script03_prefix = 'Script03_DCE_DSC/';
    script04_prefix = 'Script04_VOI/';
    script04_prefix_results = [script04_prefix 'Results/' ptno '/' ptno '_' site '/'];
    if ~isdir(script04_prefix_results)
        mkdir(script04_prefix_results)
    end

    UID = temp_series.SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(temp_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    Descrip_UID = sprintf([Descrip '_' UID]); % Descrip, SeriesInstanceUID

    if rez == 1
        img = load_nii([script01_prefix '10_Warped/Warped_' ptno '_' site '_' Descrip_UID '.nii.gz']); 
    elseif rez == 2
        img = load_nii([script01_prefix '12_Resample/Resample_' ptno '_' site '_' Descrip_UID '.nii.gz']); 
    end

    % Cylinder/Forceps
    voi_values1 = img.img(voi_idx1);
    voi_mean1 = mean(voi_values1,1);
    voi_std1 = std(voi_values1,0,1);
    voi_max1 = max(voi_values1,1);
    voi_min1 = min(voi_values1,1);
    voi_range1 = voi_max1 - voi_min1;
    voi_pix1 = numel(voi_values1);
    voi_cell1(zz,:) = {ptno site Descrip voi_mean1 voi_std1 voi_max1 voi_min1 voi_range1 voi_pix1};
    h = figure;
    hist(voi_values1)
    title([ ptno ' ' site ', ' temp_series.SeriesDescription ', mean = ' num2str(voi_mean1) ' +/- ' num2str(voi_std1)])
    saveas(h, [script04_prefix_results ptno '_' site '_Histogram_' Descrip ],'png')

    % Sphere
    voi_values2 = img.img(voi_idx2);
    voi_mean2 = mean(voi_values2,1);
    voi_std2 = std(voi_values2,0,1);
    voi_max2 = max(voi_values2,1);
    voi_min2 = min(voi_values2,1);
    voi_range2 = voi_max2 - voi_min2;
    voi_pix2 = numel(voi_values2);
    voi_cell2(zz,:) = {ptno site Descrip voi_mean2 voi_std2 voi_max2 voi_min2 voi_range2 voi_pix2};
    h = figure;
    hist(voi_values2)
    title([ ptno ' ' site ', ' temp_series.SeriesDescription ', mean = ' num2str(voi_mean2) ' +/- ' num2str(voi_std2)])
    saveas(h, [script04_prefix_results ptno '_' site '_Histogram_' Descrip ],'png')

    close all


end
