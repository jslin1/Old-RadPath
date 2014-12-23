
function HistMatch(prefix, makefile, nii_cume, n_cores, series_sql, n_intensity_levels, n_landmarks)

% Perform Histogram Matching

fid = fopen([prefix makefile],'w');
fprintf(fid, 'all:');
for ii=1:size(series_sql,1)
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:size(series_sql,1) 
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(ii));
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        prefix '10_HistMatch/Brain_HistMatch_' ptno_Descrip_UID '.nii.gz '...
        'HistogramMatch '...
        prefix '09_Truncated/Brain_Truncated0_' ptno_Descrip_UID '.nii.gz '...
        prefix nii_cume ' ' num2str(n_intensity_levels) ' ' num2str(n_landmarks)]);
end


fclose(fid);
system(['make -j ' num2str(n_cores) ' -f ' prefix makefile])







end

