

%% 0. Fix headers of DCE/DSC maps - spatial location info = applying registration works
for ii=1:size(series_DCE,1)
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_DCE(ii));
    DCE_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz']); 
    [field study]= find(cell2mat(cellfun(@(x) isequal(x,series_DCE(ii).StudyInstanceUID),struct2cell(series_DCE_maps),'UniformOutput',false)));
    for jj = study'
        [map_ptno map_Descrip map_UID map_ptno_Descrip_UID] = Generate_Label(series_DCE_maps(jj)); 
        map_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' map_ptno_Descrip_UID '.nii.gz']); 
        map_nii.hdr.dime.pixdim = DCE_nii.hdr.dime.pixdim;
        map_nii.hdr.hist = DCE_nii.hdr.hist;
        save_nii(map_nii, [map_nii.fileprefix '.nii.gz'])
    end
end

for ii=1:size(series_DSC,1)
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_DSC(ii));
    DSC_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz']); 
    [field study]= find(cell2mat(cellfun(@(x) isequal(x,series_DSC(ii).StudyInstanceUID),struct2cell(series_DSC_maps),'UniformOutput',false)));
    for jj = study'
        [map_ptno map_Descrip map_UID map_ptno_Descrip_UID] = Generate_Label(series_DSC_maps(jj)); 
        map_nii = load_nii([script03_prefix '00_radpath_raw/radpath_raw_' map_ptno_Descrip_UID '.nii.gz']); 
        map_nii.hdr.dime.pixdim = DSC_nii.hdr.dime.pixdim;
        map_nii.hdr.hist = DSC_nii.hdr.hist;
        save_nii(map_nii, [map_nii.fileprefix '.nii.gz'])
    end
end
