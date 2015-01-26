

%% 000. Copy/Re-number DICOMs for Nifti Conversion
series_sql = series_DCE_DSC_plusmaps;
for jj=1:size(series_sql,1) % for NICE/NNL, Olea DROs
    jj
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(jj));
    images   = mksqlite(['select * from Images where SeriesInstanceUID = ''' UID ''' order by SOPInstanceUID ASC']);
    pieces = strsplit(images(1).Filename,'/'); % Only works in Matlab 2014a
    halves = strsplit(images(1).Filename,pieces(end)); % Use end fragment to figure out pathname
    pathname1 = halves{1,1}; % where it's stored
    pathname2 = [script03_prefix '00_radpath_raw/DICOMs_Re-numbered/' ptno '/' ptno '_' Descrip '_DICOMs' ]; % where it's going
    if ~isdir([pathname2])
        mkdir([pathname2]) % Create new directory (if it doesn't already exist)
    end
    
    % Generate list of DICOM filenames 
    if size(dir([pathname1 '*.dcm']),1)>0
        tmp_pathname = dir([pathname1 '*.dcm']);
    elseif size(dir([pathname1 '*IM*']),1)>0
        tmp_pathname = dir([pathname1 '*IM*']);
    end
    file_names = strvcat(tmp_pathname(1:size(tmp_pathname,1)).name); % Generates file names for every file in that directory
    nimages = size(file_names,1); % # of files in that directory
    
    % Copy DICOMs and re-name based on Instance Number 
    if ~isempty(strfind(Descrip,'DCE')) || ~isempty(strfind(Descrip,'DSC')) % if it's DCE or DSC, extract last timept, and re-number
        for ii = nimages-200:nimages; 
            tmp2=dicominfo(images(ii).Filename,'dictionary','gems-dicom-dict.txt');
            if tmp2.NumberOfTemporalPositions==tmp2.TemporalPositionIdentifier % if last timept (last timept == timept of current image)
                status = copyfile(images(ii).Filename, [pathname2 '/' sprintf('%05d',tmp2.InstanceNumber)],'f'); 
            end
        end
    else % if it's a map, then copy file, and re-name based on Instance Number
        for ii = 1:nimages; 
            tmp2=dicominfo(images(ii).Filename,'dictionary','gems-dicom-dict.txt');
            status = copyfile(images(ii).Filename, [pathname2 '/' sprintf('%05d',tmp2.InstanceNumber)],'f'); % finish this command
        end
    end
end
