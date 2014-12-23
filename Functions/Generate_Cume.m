function [nii_cume] = Generate_Cume(prefix, series_sql)


for ii=1:size(series_sql,1)
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(ii));
    
    % 3a. Create cumulative histogram
    if ii==1 % For 1st image, initialize the XX_cume variable
        nii_cume = load_nii([prefix '09_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        nii_cume.img = reshape(nii_cume.img,1,[]);
    else
        nii_temp = load_nii([prefix '09_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz']);
        nii_cume.img = [nii_cume.img reshape(nii_temp.img,1,[])];
    end
    
end

imgsize_factors = factor(numel(nii_cume.img));
n3 = numel(imgsize_factors);
n1 = round(n3*1/3);
n2 = round(n3*2/3);
dim1 = prod(imgsize_factors(1:n1));
dim2 = prod(imgsize_factors(n1+1:n2));
dim3 = prod(imgsize_factors(n2+1:n3));

nii_cume.img = reshape(nii_cume.img,dim1,dim2,dim3);
nii_cume.hdr.dime.dim = [3 dim1 dim2 dim3 1 1 1 1];



end