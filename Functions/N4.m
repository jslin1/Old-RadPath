function N4(series_sql, prefix, makefile, n_cores)

    path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
    path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
    path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
    script01_prefix = 'Script01_T1_T2_SWAN/';
    script02_prefix = 'Script02_DWI_DTI/';
    script03_prefix = 'Script03_DCE_DSC/';
    script04_prefix = 'Script04_VOI/';
    LPBA40_Template_location= '/mnt/data/scratch/igilab/jslin1/LPBA40_Template.nii.gz';
    LPBA40_mask_location = '/mnt/data/scratch/igilab/jslin1/LPBA40_mask.nii.gz';
    mksqlite('open', 'ctkDICOM.sql' );

fid = fopen([prefix  makefile],'w');
fprintf(fid, 'all:');
for ii=1:size(series_sql,1)
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:size(series_sql,1)
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(ii));

    % 1. N4 Bias Field Correction
    fprintf(fid,['\t$(ANTSPATH)/N4BiasFieldCorrection -d 3 '...
		'-i ' script01_prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz -s 4 -c [50x50x50x50,1e-7] -b [200] '...
		'-o ' script01_prefix '01_N4/N4_' ptno_Descrip_UID '.nii.gz\n']);
end

fclose(fid);
system(['make -j ' num2str(n_cores) ' -f ' prefix makefile])


end


