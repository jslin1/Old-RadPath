function BO(lot, prefix, mask_series, input_series)

path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql

for ii=1:size(mask_series,1)

[m_ptno m_Descrip m_UID m_ptno_Descrip_UID] = Generate_Label(mask_series(ii));
[i_field i_study]= find(cell2mat(cellfun(@(x) isequal(x,mask_series(ii).StudyInstanceUID),struct2cell(input_series),'UniformOutput',false)));

    for jj = i_study'
        [i_ptno i_Descrip i_UID i_ptno_Descrip_UID] = Generate_Label(input_series(jj));
        
        if lot==0
            m_dir = [script01_prefix '04_Masks/FINALMASK_'];
            i_dir = [prefix '01_N4/N4_'];
            o_dir = [prefix '08_Brain_Only/Brain_'];
        elseif lot==1 || lot==2
            m_dir = [script01_prefix '04_Masks/FINALMASK_'];
            i_dir = [prefix '07_Warped/Warped_'];
            o_dir = [prefix '08_Brain_Only/Brain_'];
        elseif lot==3
            m_dir = [prefix '04_Masks/FINALMASK_'];
            if ~isempty(strfind(i_ptno_Descrip_UID,'DCE')) || ~isempty(strfind(i_ptno_Descrip_UID,'DSC'))
                i_dir = [prefix '01_N4/N4_'];
            else
                i_dir = [prefix '00_radpath_raw/radpath_raw_'];
            end
            o_dir = [prefix '05_Brain_Only/Brain_']; 
        end

        m_title = [m_dir m_ptno_Descrip_UID '.nii.Labelfield.nii'];
        i_title = [i_dir i_ptno_Descrip_UID '.nii.gz'];
        o_title = [o_dir i_ptno_Descrip_UID '.nii.gz'];
 
        m_file = load_nii(m_title);
        i_file = load_nii(i_title);
        i_file.img = double(i_file.img).*double(m_file.img);
        save_nii(i_file, o_title)
        
        
% o_file = load_nii(o_title);
% figure
% imagesc(e_file.img(:,:,8))
% figure
% imagesc(o_file.img(:,:,8))
        
% Flips y axis of pictures so not using this command
%         system(['/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/MultiplyImages 3 '...
%             m_title ' '...
%             e_title ' '...
%             o_title '']);
        
    end
end

end
