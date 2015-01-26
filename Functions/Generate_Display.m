disp_side = 100; % side length in mm of extracted portion of image

% World --> Vox
% World 1 --> Vox 1/2/3
if srow_x(1)~=0 
    i = 1 + (-point(1) - srow_x(4))/srow_x(1); % ITK i
    origin_i = double(round(i-(disp_side/abs(srow_x(1)))/2));
    num_i = double(round(disp_side/abs(srow_x(1))));
elseif srow_x(2)~=0 
    j = 1 + (-point(1) - srow_x(4))/srow_x(2); % ITK j
    origin_j = double(round(j-(disp_side/abs(srow_x(2)))/2));
    num_j = double(round(disp_side/abs(srow_x(2))));
elseif srow_x(3)~=0 
    k = 1 + (-point(1) - srow_x(4))/srow_x(3); % ITK k
    origin_k = double(round(k-(disp_side/abs(srow_x(3)))/2));
    num_k = double(round(disp_side/abs(srow_x(3))));
end

% World 2, Vox 1/2/3
if srow_y(1)~=0 
    i = 1 + (-point(2) - srow_y(4))/srow_y(1); % ITK i
    origin_i = double(round(i-(disp_side/abs(srow_y(1)))/2));
    num_i = double(round(disp_side/abs(srow_y(1))));
elseif srow_y(2)~=0 
    j = 1 + (-point(2) - srow_y(4))/srow_y(2); % ITK i
    origin_j = double(round(j-(disp_side/abs(srow_y(2)))/2));
    num_j = double(round(disp_side/abs(srow_y(2))));
elseif srow_y(3)~=0 
    k = 1 + (-point(2) - srow_y(4))/srow_y(3); % ITK i
    origin_k = double(round(k-(disp_side/abs(srow_y(3)))/2));
    num_k = double(round(disp_side/abs(srow_y(3))));
end

% World 3, Vox 1/2/3
if srow_z(1)~=0 
    i = 1 + (point(3) - srow_z(4))/srow_z(1); % ITK j
    origin_i = double(round(i-(disp_side/abs(srow_z(1)))/2));
    num_i = double(round(disp_side/abs(srow_z(1))));
elseif srow_z(2)~=0 
    j = 1 + (point(3) - srow_z(4))/srow_z(2); % ITK j
    origin_j = double(round(j-(disp_side/abs(srow_z(2)))/2));
    num_j = double(round(disp_side/abs(srow_z(2))));
elseif srow_z(3)~=0 
    k = 1 + (point(3) - srow_z(4))/srow_z(3); % ITK j
    origin_k = double(round(k-(disp_side/abs(srow_z(3)))/2));
    num_k = double(round(disp_side/abs(srow_z(3))));
end 

 extract_origin = [num2str(origin_i) 'x' num2str(origin_j) 'x' num2str(origin_k)];
 extract_size = [num2str(num_i) 'x' num2str(num_j) 'x' num2str(num_k)];

% Series for which I want to generate display volumes
series_list = {series_T2_register      series_ADC      series_rBV      series_Ktrans1   series_Vp1      series_AUC1};
for kk = 1:size(series_list,2)
    d_series = series_list{:,kk};
    [d_field d_study]= find(cell2mat(cellfun(@(x) isequal(x,studies_all(jj).StudyInstanceUID),struct2cell(d_series),'UniformOutput',false)));
    [d_ptno d_Descrip d_UID d_ptno_Descrip_UID] = Generate_Label(d_series(d_study));
    
    if ~isempty(strfind(d_Descrip,'T1')) || ~isempty(strfind(d_Descrip,'T2')) || ~isempty(strfind(d_Descrip,'SWAN'))
        prefix = script01_prefix;
        d_conv_factor = double(1);
        source_dir = ['07_Warped/Warped_'];
%         source_dir = ['09_GaussNorm/GN_'];
        if kk==1 % T2 = always the 1st one
            source_dir = ['01_N4/N4_'];
        end
    elseif ~isempty(strfind(d_Descrip,'Apparent_Diffusion_Coefficient')) || ~isempty(strfind(d_Descrip,'Fractional_Aniso')) || ~isempty(strfind(d_Descrip,'Average_DC'))
        prefix = script02_prefix;
        d_conv_factor = double(1e-6);
        source_dir = ['07_Warped/Warped_'];
%         source_dir = ['09_GaussNorm/GN_'];
    elseif ~isempty(strfind(d_Descrip,'DSC')) || ~isempty(strfind(d_Descrip,'DCE'))...
            || ~isempty(strfind(d_Descrip,'rBV')) || ~isempty(strfind(d_Descrip,'rBF'))...
            || ~isempty(strfind(d_Descrip,'MTT')) || ~isempty(strfind(d_Descrip,'Delay'))...
            || ~isempty(strfind(d_Descrip,'Leakage__K2_')) || ~isempty(strfind(d_Descrip,'K12'))...
            || ~isempty(strfind(d_Descrip,'K21')) || ~isempty(strfind(d_Descrip,'Vp'))...
            || ~isempty(strfind(d_Descrip,'Ve__distribution volume_')) || ~isempty(strfind(d_Descrip,'Wash-in'))...
            || ~isempty(strfind(d_Descrip,'Wash-out')) || ~isempty(strfind(d_Descrip,'Time_to_peak'))...
            || ~isempty(strfind(d_Descrip,'Area_under_curve')) || ~isempty(strfind(d_Descrip,'Peak_enhancement'))
        prefix = script03_prefix;
        images = mksqlite(['select * from Images where SeriesInstanceUID = ''' d_UID ''' order by Filename ASC' ]);
        tmp2=dicominfo(images(1).Filename,'dictionary','gems-dicom-dict.txt');
        d_conv_factor = double(tmp2.Private_0077_1001);
        source_dir = ['08_Warped/Warped_'];
%         source_dir = ['09_GaussNorm/GN_'];
    end

    system(['c3d ' prefix source_dir d_ptno_Descrip_UID '.nii.gz '... % Input filename
            '-region ' extract_origin 'vox ' extract_size 'vox '... % Origin-ijk dimensions-ijk
             '-type double -noround '...
            '-o ' prefix '10_Extract/Extract_' ptno_site '_' d_ptno_Descrip_UID '.nii.gz']);

    system(['c3d ' prefix '10_Extract/Extract_' ptno_site '_' d_ptno_Descrip_UID '.nii.gz '... % Input filename
            '-interpolation NearestNeighbor -resample-mm 0.2x0.2x0.2mm '... % Interpolation type, resampling dimensions
            '-type double -noround '...
            '-o ' prefix '11_Resample/Resample_' ptno_site '_' d_ptno_Descrip_UID '.nii.gz']);
        
    d_temp = load_nii([prefix source_dir d_ptno_Descrip_UID '.nii.gz']);
    d_values = double(d_temp.img(d_temp.img>0))*d_conv_factor;
    [d_mean d_std d_min d_max d_range d_numel]=Generate_Stats(d_values);
    
    figure;		hist(d_values);        xlabel('Values');        ylabel('# Voxels')
    title({     [ d_ptno ', ' d_Descrip];...
                ['mean = '      num2str(d_mean)...
                ' +/- '         num2str(d_std)...
                ', # vox = '    num2str(d_numel)]}, 'Interpreter' ,'none')
    saveas(gcf, [script04_results ptno '/' 'Histogram_' d_ptno '_' Descrip],'png')
    close all

end
    
    
    