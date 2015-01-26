function [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(temp_series)

path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql

UID = temp_series.SeriesInstanceUID;
Descrip = strrep(strrep(strrep(strrep(strrep(strrep(temp_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star'),'.','');
images = mksqlite(['select * from Images where SeriesInstanceUID = ''' UID ''' order by Filename ASC' ]);
tmp2=dicominfo(images(1).Filename,'dictionary','gems-dicom-dict.txt');
ptno = mrn_to_ptno(tmp2.PatientID);
ptno_Descrip_UID = sprintf([ptno '_' Descrip '_' UID]); % patient #, Descrip, SeriesInstanceUID

end