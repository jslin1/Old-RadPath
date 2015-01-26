function Create_Mask(lot, prefix, makefile, n_cores, series_sql, MD_parameter, MC_parameter)

path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql
LPBA40_Template_location= '/mnt/data/scratch/igilab/jslin1/LPBA40_Template.nii.gz';
LPBA40_mask_location = '/mnt/data/scratch/igilab/jslin1/LPBA40_mask.nii.gz';

if lot==0 || lot==1 || lot==3
    title_pre = [prefix '01_N4/N4_' ];
elseif lot==2
    title_pre = [prefix '00_radpath_raw/radpath_raw_' ];
end


fid = fopen([prefix  makefile],'w');
fprintf(fid, 'all:');
for ii=1:size(series_sql,1)
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:size(series_sql,1)
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    [ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(ii));
    
    title = [title_pre ptno_Descrip_UID '.nii.gz'];
    m_rinput_dir = [prefix '02_Mask_RInput/'];
    m_routput_dir = [prefix '03_Mask_ROutput/'];
    m_dir = [prefix '04_Masks/'];
    
    % 1. Perform rough down-sampled registration as beginning step to help initialize later registration
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
 		LPBA40_Template_location ' '...
        m_rinput_dir 'Affine_Fixed_LPBA40.nii.gz 4 4 4 1\n']);
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
		title ' '...
		m_rinput_dir 'Affine_Moving_' ptno_Descrip_UID '.nii.gz 4 4 4 1\n']);
    fprintf(fid,['\t$(ANTSPATH)/antsAffineInitializer 3 '...
		m_rinput_dir 'Affine_Fixed_LPBA40.nii.gz '...
		m_rinput_dir 'Affine_Moving_' ptno_Descrip_UID '.nii.gz '...
	    m_rinput_dir 'InitialAffine_' ptno_Descrip_UID '.mat 15 0.1 0 10\n']);

    % 2. Create Laplacian (edge detection) version of images to assist in registration
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		m_rinput_dir 'Laplacian_' ptno_Descrip_UID '.nii.gz '...
		'Laplacian '...
		title ' 1.5 1\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		m_rinput_dir 'Laplacian_LPBA40.nii.gz '...
		'Laplacian '...
		LPBA40_Template_location ' 1.5 1\n']);

    % 3. Perform registration of image to template
    fprintf(fid,['\t$(ANTSPATH)/antsRegistration -d 3 -u 1 -w [0.025, 0.975] '...
	   '-o ' m_routput_dir 'RegistrationOutput_' ptno_Descrip_UID '_ '...
	   '-r ' m_rinput_dir 'InitialAffine_' ptno_Descrip_UID '.mat -z 1 --float 0 '...
        '-m MI[' LPBA40_Template_location ',' title ',1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m MI[' LPBA40_Template_location ',' title ',1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m CC[' LPBA40_Template_location ',' title ',0.5,4] '...
        '-m CC[' m_rinput_dir 'Laplacian_LPBA40.nii.gz,'...
                 m_rinput_dir 'Laplacian_' ptno_Descrip_UID '.nii.gz,0.5,4] -c [50x10x0,1e-9,15] -t SyN[0.1,3,0] -f 4x2x1 -s 2x1x0\n']);

    % 4. Apply registration to mask
    fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
		'-i ' LPBA40_mask_location ' '...
		'-o ' m_dir 'Mask_' ptno_Descrip_UID '.nii.gz '...
		'-r ' title ' -n Gaussian '...
        '-t [' m_routput_dir 'RegistrationOutput_' ptno_Descrip_UID '_0GenericAffine.mat,1] '...
		'-t '  m_routput_dir 'RegistrationOutput_' ptno_Descrip_UID '_1InverseWarp.nii.gz --float 0\n']);

    % 5. Treshold and Get largest component of mask. 
    fprintf(fid,['\t$(ANTSPATH)/ThresholdImage 3 '...
		m_dir 'Mask_' ptno_Descrip_UID '.nii.gz '...
		m_dir 'Mask_thresholded_' ptno_Descrip_UID '.nii.gz 0.5 1 1 0\n']); % Threshold Mask
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		m_dir 'Mask_LC_' ptno_Descrip_UID '.nii.gz '...
		'GetLargestComponent '...
		m_dir 'Mask_thresholded_' ptno_Descrip_UID '.nii.gz\n']); % Get largest component of mask 

    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		m_dir 'Mask_MD_' ptno_Descrip_UID '.nii.gz '...
		'MD '...
		m_dir 'Mask_LC_' ptno_Descrip_UID '.nii.gz ' num2str(MD_parameter) '\n']); % Dilate mask
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		m_dir 'Mask_MC_' ptno_Descrip_UID '.nii.gz '...
		'MC '...
		m_dir 'Mask_MD_' ptno_Descrip_UID '.nii.gz ' num2str(MC_parameter) '\n']); % Close mask
    
    % 6. Perform Atropos segmentation to get WM, GM, and CSF. Combine these masks.
    fprintf(fid,['\t$(ANTSPATH)/Atropos -d 3 '...
		'-o ' m_dir 'Atropos_Mask_' ptno_Descrip_UID '.nii.gz '...
		'-a ' title ' '...
		'-x ' m_dir 'Mask_MC_' ptno_Descrip_UID '.nii.gz '...
        	'-i kmeans[3] -c [3, 0.0] -m [0.1, 1x1x1] -k Gaussian\n']);

    % 7. Separate WM, GM, CSF masks
    fprintf(fid,['\t$(ANTSPATH)/ThresholdImage 3 '...
		m_dir 'Atropos_Mask_' ptno_Descrip_UID '.nii.gz '...
		m_dir 'Atropos_WM_' ptno_Descrip_UID '.nii.gz 3 3 1 0\n']);
    fprintf(fid,['\t$(ANTSPATH)/ThresholdImage 3 '...
		m_dir 'Atropos_Mask_' ptno_Descrip_UID '.nii.gz '...
		m_dir 'Atropos_GM_' ptno_Descrip_UID '.nii.gz 2 2 1 0\n']);
    fprintf(fid,['\t$(ANTSPATH)/ThresholdImage 3 '...
		m_dir 'Atropos_Mask_' ptno_Descrip_UID '.nii.gz '...
		m_dir 'Atropos_CSF_' ptno_Descrip_UID '.nii.gz 1 1 1 0\n']);

    % 8. Combine WM, GM, CSF masks into 1
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		m_dir 'Atropos_WM_GM_' ptno_Descrip_UID '.nii.gz '...
		'addtozero '...
		m_dir 'Atropos_WM_' ptno_Descrip_UID '.nii.gz '...
		m_dir 'Atropos_GM_' ptno_Descrip_UID '.nii.gz\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		m_dir 'Atropos_WM_GM_CSF_' ptno_Descrip_UID '.nii.gz '...
		'addtozero '...
		m_dir 'Atropos_WM_GM_' ptno_Descrip_UID '.nii.gz '...
		m_dir 'Atropos_CSF_' ptno_Descrip_UID '.nii.gz\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		m_dir 'Atropos_LC_' ptno_Descrip_UID '.nii.gz '...
		'GetLargestComponent '...
		m_dir 'Atropos_WM_GM_CSF_' ptno_Descrip_UID '.nii.gz\n']);

end
fclose(fid);
system(['make -j ' num2str(n_cores) ' -f ' prefix makefile])

end


