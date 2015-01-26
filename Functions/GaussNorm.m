function GaussNorm(lot, prefix, series_sql)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if lot==0 || lot==1 || lot==2
    BO_dir = [prefix '08_Brain_Only/Brain_'];
    GN_dir = [prefix '09_GaussNorm/GN_'];
elseif lot==3
    BO_dir = [prefix '08_Warped/Warped_'];
    GN_dir = [prefix '09_GaussNorm/GN_'];
end
    
for ii=1:size(series_sql,1)
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(ii));
    i_title = [BO_dir ptno_Descrip_UID '.nii.gz'];
    o_title = [GN_dir ptno_Descrip_UID '.nii.gz'];

    file = load_nii(i_title);
    brain_values = double(file.img(file.img(:)>0));
    brain_std = std(brain_values,0,1);
    file.img = double(file.img/brain_std);
    save_nii(file, o_title)
end

end

