function DICOM_to_Nifti(lot, prefix, makefile, series_sql, n_cores)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fid = fopen([prefix makefile],'w');
fprintf(fid, 'all:');
for ii=1:size(series_sql,1)
    fprintf(fid, [' job' num2str(ii)]);
end

for jj = 1:size(series_sql,1) 
    fprintf(fid,['\n\njob' num2str(jj) ':\n']);
    
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(jj));
    if lot==0 || lot==1 || lot==2
        images   = mksqlite(['select * from Images where SeriesInstanceUID = ''' UID '''' ]);
        pieces = strsplit(images(1).Filename,'/'); % Only works in Matlab 2014a
        halves = strsplit(images(1).Filename,pieces(end)); % Use end fragment to figure out pathname
        pathname = strrep(halves{1,1},' ','_');
    elseif lot==3
        pathname = [prefix '00_radpath_raw/DICOMs_Re-numbered/' ptno '/' ptno '_' Descrip '_DICOMs' ];    
    end
    o_title = [prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz'];
    
    % Convert 1 series to Nifti - Command Source_Directory Output_Filename SeriesInstanceUID(to ensure only that series is used)
    fprintf(fid, ['\tDicomSeriesReadImageWrite2 ' ...
        pathname ' '...
        o_title ' '...
        UID '\n']);
end

fclose(fid);
system(['make -j ' num2str(n_cores) ' -f ' prefix makefile])

end

