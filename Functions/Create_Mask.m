function Create_Mask(fid, ii, series_sql_list, MD_parameter, MC_parameter)

    path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
    path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
    path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
    script01_prefix = 'Script01_T1_T2_SWAN/';
    script02_prefix = 'Script02_DWI_DTI/';
    script03_prefix = 'Script03_DCE_DSC/';
    script04_prefix = 'Script04_VOI/';
    LPBA40_Template_location= '/mnt/data/scratch/igilab/jslin1/LPBA40_Template.nii.gz';
    LPBA40_mask_location = '/mnt/data/scratch/igilab/jslin1/LPBA40_mask.nii.gz';

    UID = series_sql_list(ii).SeriesInstanceUID;
    Descrip = strrep(strrep(strrep(strrep(strrep(series_sql_list(ii).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    ptno_Descrip_UID = sprintf(['%02d_' Descrip '_' UID],ii); %

%%%%%%%%%%%%%%%%%%%%%% N4, Brain Extract, Truncation, HistMatch
    % 1. N4 Bias Field Correction
    fprintf(fid,['\t$(ANTSPATH)/N4BiasFieldCorrection -d 3 '...
		'-i ' script01_prefix '00_radpath_raw/radpath_raw_' ptno_Descrip_UID '.nii.gz -s 4 -c [50x50x50x50,1e-7] -b [200] '...
		'-o ' script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz\n']);

    % 2. Create Mask and Extract Brain
    % 2a. Perform rough down-sampled registration as beginning step to help initialize later registration
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
 		LPBA40_Template_location ' '...
		script01_prefix '02_RInput/Affine_Fixed_LPBA40.nii.gz 4 4 4 1\n']);
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '02_RInput/Affine_Moving_' ptno_Descrip_UID '.nii.gz 4 4 4 1\n']);
    fprintf(fid,['\t$(ANTSPATH)/antsAffineInitializer 3 '...
		script01_prefix '02_RInput/Affine_Fixed_LPBA40.nii.gz '...
		script01_prefix '02_RInput/Affine_Moving_' ptno_Descrip_UID '.nii.gz '...
	    script01_prefix '02_RInput/InitialAffine_' ptno_Descrip_UID '.mat 15 0.1 0 10\n']);

    % 2b. Create Laplacian (edge detection) version of images to assist in registration
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '02_RInput/Laplacian_' ptno_Descrip_UID '.nii.gz '...
		'Laplacian '...
		script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz 1.5 1\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '02_RInput/Laplacian_LPBA40.nii.gz '...
		'Laplacian '...
		 LPBA40_Template_location ' 1.5 1\n']);

    % 2c. Perform registration of image to template
    fprintf(fid,['\t$(ANTSPATH)/antsRegistration -d 3 -u 1 -w [0.025, 0.975] '...
	   '-o ' script01_prefix '03_ROutput/RegistrationOutput_' ptno_Descrip_UID '_ '...
	   '-r ' script01_prefix '02_RInput/InitialAffine_' ptno_Descrip_UID '.mat -z 1 --float 0 '...
        '-m MI[' LPBA40_Template_location ','...
	      script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m MI[' LPBA40_Template_location ','...
	      script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m CC[' LPBA40_Template_location ','...
	      script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz,0.5,4] '...
        '-m CC[02_RInput/Laplacian_LPBA40.nii.gz,'...
	      script01_prefix '02_RInput/Laplacian_' ptno_Descrip_UID '.nii.gz,0.5,4] -c [50x10x0,1e-9,15] -t SyN[0.1,3,0] -f 4x2x1 -s 2x1x0\n']);

    % 2d. Apply registration to mask
    fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
		'-i ' LPBA40_mask_location ' '...
		'-o ' script01_prefix '04_Masks/Mask_' ptno_Descrip_UID '.nii.gz '...
		'-r ' script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz -n Gaussian '...
        '-t [' script01_prefix '03_ROutput/RegistrationOutput_' ptno_Descrip_UID '_0GenericAffine.mat,1] '...
		'-t ' script01_prefix '03_ROutput/RegistrationOutput_' ptno_Descrip_UID '_1InverseWarp.nii.gz --float 0\n']);

    % 2e. Treshold and Get largest component of mask. 
    fprintf(fid,['\t$(ANTSPATH)/ThresholdImage 3 '...
		script01_prefix '04_Masks/Mask_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '04_Masks/Mask_thresholded_' ptno_Descrip_UID '.nii.gz 0.5 1 1 0\n']); % Threshold Mask
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '04_Masks/Mask_LC_' ptno_Descrip_UID '.nii.gz '...
		'GetLargestComponent '...
		script01_prefix '04_Masks/Mask_thresholded_' ptno_Descrip_UID '.nii.gz\n']); % Get largest component of mask 

    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '04_Masks/Mask_MD_' ptno_Descrip_UID '.nii.gz '...
		'MD '...
		script01_prefix '04_Masks/Mask_LC_' ptno_Descrip_UID '.nii.gz ' num2str(MD_parameter) '\n']); % Dilate mask
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '04_Masks/Mask_MC_' ptno_Descrip_UID '.nii.gz '...
		'MC '...
		script01_prefix '04_Masks/Mask_MD_' ptno_Descrip_UID '.nii.gz ' num2str(MC_parameter) '\n']); % Close mask
    
    % 2f. Perform Atropos segmentation to get WM, GM, and CSF. Combine these masks.
    fprintf(fid,['\t$(ANTSPATH)/Atropos -d 3 '...
		'-o ' script01_prefix '04_Masks/Atropos_Mask_' ptno_Descrip_UID '.nii.gz '...
		'-a ' script01_prefix '01_N4/Brain_N4_' ptno_Descrip_UID '.nii.gz '...
		'-x ' script01_prefix '04_Masks/Mask_MC_' ptno_Descrip_UID '.nii.gz '...
        	'-i kmeans[3] -c [3, 0.0] -m [0.1, 1x1x1] -k Gaussian\n']);

    % 2g. Separate WM, GM, CSF masks
    fprintf(fid,['\t$(ANTSPATH)/ThresholdImage 3 '...
		script01_prefix '04_Masks/Atropos_Mask_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '04_Masks/Atropos_WM_' ptno_Descrip_UID '.nii.gz 3 3 1 0\n']);
    fprintf(fid,['\t$(ANTSPATH)/ThresholdImage 3 '...
		script01_prefix '04_Masks/Atropos_Mask_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '04_Masks/Atropos_GM_' ptno_Descrip_UID '.nii.gz 2 2 1 0\n']);
    fprintf(fid,['\t$(ANTSPATH)/ThresholdImage 3 '...
		script01_prefix '04_Masks/Atropos_Mask_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '04_Masks/Atropos_CSF_' ptno_Descrip_UID '.nii.gz 1 1 1 0\n']);

    % 2h. Combine WM, GM, CSF masks into 1
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '04_Masks/Atropos_WM_GM_' ptno_Descrip_UID '.nii.gz '...
		'addtozero '...
		script01_prefix '04_Masks/Atropos_WM_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '04_Masks/Atropos_GM_' ptno_Descrip_UID '.nii.gz\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '04_Masks/Atropos_WM_GM_CSF_' ptno_Descrip_UID '.nii.gz '...
		'addtozero '...
		script01_prefix '04_Masks/Atropos_WM_GM_' ptno_Descrip_UID '.nii.gz '...
		script01_prefix '04_Masks/Atropos_CSF_' ptno_Descrip_UID '.nii.gz\n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
		script01_prefix '04_Masks/Atropos_LC_' ptno_Descrip_UID '.nii.gz '...
		'GetLargestComponent '...
		script01_prefix '04_Masks/Atropos_WM_GM_CSF_' ptno_Descrip_UID '.nii.gz\n']);

end


