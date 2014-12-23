function Register3(prefix, makefile, n_cores, fixed_series, moving_series, apply_series)

path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';
mksqlite('open', 'ctkDICOM.sql');

fid = fopen([prefix makefile],'w');
fprintf(fid, 'all:');
for ii=1:size(moving_series,1)
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:size(moving_series,1)
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);

    [Moving_ptno Moving_Descrip Moving_UID Moving_ptno_Descrip_UID] = Generate_Label(moving_series(ii));
    [fixed_field fixed_study]= find(cell2mat(cellfun(@(x) isequal(x,moving_series(ii).StudyInstanceUID),struct2cell(fixed_series),'UniformOutput',false)));
    [Fixed_ptno Fixed_Descrip Fixed_UID Fixed_ptno_Descrip_UID] = Generate_Label(fixed_series(fixed_study));

% Add lasttimept as text
% Use AUC?

    % 1. Perform rough down-sampled registration as beginning step to help
    % initialize later registration
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        script01_prefix '01_N4/N4_' Fixed_ptno_Descrip_UID '.nii.gz '...
        prefix '02_RInput/Affine_Fixed_' Fixed_ptno_Descrip_UID '.nii.gz 4 4 4 1\n']); % Down-sample fixed image
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        prefix '00_radpath_raw/radpath_raw_' Moving_ptno_Descrip_UID '.nii.gz '...
        prefix '02_RInput/Affine_Moving_' Moving_ptno_Descrip_UID '.nii.gz 4 4 4 1\n']); %Down-sample moving image
    fprintf(fid,['\t$(ANTSPATH)/antsAffineInitializer 3 '...
        prefix '02_RInput/Affine_Fixed_' Fixed_ptno_Descrip_UID '.nii.gz '...
        prefix '02_RInput/Affine_Moving_' Moving_ptno_Descrip_UID '.nii.gz '...
        prefix '02_RInput/InitialAffine_' Moving_ptno_Descrip_UID '.mat 15 0.1 0 10\n']); % Rough registration

    % 2. Create Laplacian (edge detection) version of images to assist in registration
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        prefix '02_RInput/Laplacian_' Fixed_ptno_Descrip_UID '.nii.gz '...
        'Laplacian '...
        script01_prefix '01_N4/N4_' Fixed_ptno_Descrip_UID '.nii.gz \n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        prefix '02_RInput/Laplacian_' Moving_ptno_Descrip_UID '.nii.gz '...
        'Laplacian '...
        prefix '00_radpath_raw/radpath_raw_' Moving_ptno_Descrip_UID '.nii.gz 1.5 1\n']);

    % 3. Perform registration 
    fprintf(fid,['\t$(ANTSPATH)/antsRegistration -d 3 -u 1 -w [0.025, 0.975] '...
        '-o ' prefix '03_ROutput/RegistrationOutput_' Moving_ptno_Descrip_UID '_ '...
        '-r ' prefix '02_RInput/InitialAffine_' Moving_ptno_Descrip_UID '.mat -z 1 --float 0 '...
        '-m MI[' script01_prefix '01_N4/N4_' Fixed_ptno_Descrip_UID '.nii.gz,'...
              prefix '00_radpath_raw/radpath_raw_' Moving_ptno_Descrip_UID '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m MI[' script01_prefix '01_N4/N4_' Fixed_ptno_Descrip_UID '.nii.gz,'...
              prefix '00_radpath_raw/radpath_raw_' Moving_ptno_Descrip_UID '.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
              '>' prefix '03_ROutput/RegistrationOutput_' Moving_ptno_Descrip_UID '.txt \n']);

    % 4. Apply registration     
    [apply_field apply_study]= find(cell2mat(cellfun(@(x) isequal(x,moving_series(ii).StudyInstanceUID),struct2cell(apply_series),'UniformOutput',false)));
    for jj = apply_study
        [Apply_ptno Apply_Descrip Apply_UID Apply_ptno_Descrip_UID] = Generate_Label(apply_series(jj));

        fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
                    '-i ' prefix '00_radpath_raw/radpath_raw_' Apply_ptno_Descrip_UID '.nii.gz '...
                    '-o ' prefix '04_Warped/Warped_' Apply_ptno_Descrip_UID '.nii.gz '... 
                    '-r ' script01_prefix '01_N4/N4_' Fixed_ptno_Descrip_UID '.nii.gz '...
                    '-n NearestNeighbor '... % Interpolation method for resampling
                    '-t [' prefix '03_ROutput/RegistrationOutput_' Moving_ptno_Descrip_UID '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
                    '\n']); 
    end
end

fclose(fid);
system(['make -j ' num2str(n_cores) ' -f ' prefix makefile])




end