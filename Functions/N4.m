function N4(series_sql, prefix, makefile, n_cores)

path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql

fid = fopen([prefix  makefile],'w');
fprintf(fid, 'all:');
for ii=1:size(series_sql,1)
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:size(series_sql,1)
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(ii));
    i_title = [prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz'];
    o_title = [prefix '01_N4/N4_' ptno_Descrip_UID '.nii.gz'];
    
    % 1. N4 Bias Field Correction
    fprintf(fid,['\t$(ANTSPATH)/N4BiasFieldCorrection -d 3 '...
		'-i ' i_title ' -s 4 -c [50x50x50x50,1e-7] -b [200] '...
		'-o ' o_title '\n']);
end

fclose(fid);
system(['make -j ' num2str(n_cores) ' -f ' prefix makefile])


end


