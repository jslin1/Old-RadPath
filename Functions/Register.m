function Register(lot, prefix, makefile, n_cores, fixed_series, moving_series, laplace_series, apply_series)

path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql

f_dir = [script01_prefix '01_N4/N4_'];
% f_dir = [script01_prefix '08_Brain_Only/Brain_'];

if lot==1 % t1,t2,swan
    m_dir = [prefix '01_N4/N4_' ];
    rinput_dir = [prefix '05_RInput/'];
    routput_dir = [prefix '06_ROutput/'];
    w_dir = [prefix '07_Warped/Warped_'];
    a_dir = m_dir;
elseif lot==2
    m_dir = [prefix '00_radpath_raw/radpath_raw_' ];
    rinput_dir = [prefix '05_RInput/'];
    routput_dir = [prefix '06_ROutput/'];
    w_dir = [prefix '07_Warped/Warped_'];
    a_dir = m_dir;
elseif lot==3 % dsc, dce
    m_dir = [prefix '05_Brain_Only/Brain_' ];
    rinput_dir = [prefix '06_RInput/'];
    routput_dir = [prefix '07_ROutput/'];
    w_dir = [prefix '08_Warped/Warped_'];
    a_dir = m_dir;
end

fid = fopen([prefix makefile],'w');
fprintf(fid, 'all:');
for ii=1:size(moving_series,1)
    fprintf(fid, [' job' num2str(ii)]);
end
fprintf(fid, '\nANTSPATH=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin\n');

for ii=1:size(moving_series,1)
    fprintf(fid,['\n\njob' num2str(ii) ':\n']);

    [m_ptno m_Descrip m_UID m_ptno_Descrip_UID] = Generate_Label(moving_series(ii));
    [f_field f_study]= find(cell2mat(cellfun(@(x) isequal(x,moving_series(ii).StudyInstanceUID),struct2cell(fixed_series),'UniformOutput',false)));
    [f_ptno f_Descrip f_UID f_ptno_Descrip_UID] = Generate_Label(fixed_series(f_study));

    [l_field l_study]= find(cell2mat(cellfun(@(x) isequal(x,moving_series(ii).StudyInstanceUID),struct2cell(laplace_series),'UniformOutput',false)));
    [l_ptno l_Descrip l_UID l_ptno_Descrip_UID] = Generate_Label(laplace_series(l_study));
    
    
    f_title = [f_dir f_ptno_Descrip_UID '.nii.gz'];
    m_title = [m_dir m_ptno_Descrip_UID '.nii.gz'];
    l_title = [m_dir l_ptno_Descrip_UID '.nii.gz'];
    f_Laplacian = [rinput_dir 'Laplacian_' f_ptno_Descrip_UID '.nii.gz'];
    m_Laplacian = [rinput_dir 'Laplacian_' m_ptno_Descrip_UID '.nii.gz'];
    l_Laplacian = [rinput_dir 'Laplacian_' l_ptno_Descrip_UID '.nii.gz'];
    
    % 1. Perform rough down-sampled registration as beginning step to help
    % initialize later registration
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        f_title ' '...
        rinput_dir 'Affine_Fixed_' f_ptno_Descrip_UID '.nii.gz 4 4 4 1\n']); % Down-sample fixed image
    fprintf(fid,['\t$(ANTSPATH)/ResampleImageBySpacing 3 '...
        m_title ' '...
        rinput_dir 'Affine_Moving_' m_ptno_Descrip_UID '.nii.gz 4 4 4 1\n']); %Down-sample moving image
    fprintf(fid,['\t$(ANTSPATH)/antsAffineInitializer 3 '...
        rinput_dir 'Affine_Fixed_' f_ptno_Descrip_UID '.nii.gz '...
        rinput_dir 'Affine_Moving_' m_ptno_Descrip_UID '.nii.gz '...
        rinput_dir 'InitialAffine_' m_ptno_Descrip_UID '.mat 15 0.1 0 10\n']); % Rough registration

    % 2. Create Laplacian (edge detection) version of images to assist in registration
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        f_Laplacian ' '...
        'Laplacian '...
        f_title ' \n']);
    fprintf(fid,['\t$(ANTSPATH)/ImageMath 3 '...
        l_Laplacian ' '...
        'Laplacian '...
        l_title ' 1.5 1\n']);

    % 3. Perform registration 
    fprintf(fid,['\t$(ANTSPATH)/antsRegistration -d 3 -u 1 -w [0.025, 0.975] '...
        '-o ' routput_dir 'RegistrationOutput_' m_ptno_Descrip_UID '_ '...
        '-r ' rinput_dir 'InitialAffine_' m_ptno_Descrip_UID '.mat -z 1 --float 0 '...
        '-m MI[' f_title ',' m_title ',1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Rigid[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m MI[' f_title ',' m_title ',1,32,Regular,0.25] -c [1000x500x250x100,1e-8,10] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
        '-m MI[' f_Laplacian ',' l_Laplacian ',1,32,Regular,0.25] '...
                '-c [1000x500x250x100,1e-8,10] -t Affine[0.1] -f 8x4x2x1 -s 4x2x1x0 '...
                '>' routput_dir 'RegistrationOutput_' m_ptno_Descrip_UID '.txt \n']);
            
    % 4. Apply registration 
    [a_field a_study]= find(cell2mat(cellfun(@(x) isequal(x,moving_series(ii).StudyInstanceUID),struct2cell(apply_series),'UniformOutput',false)));
    for jj = a_study'
        [a_ptno a_Descrip a_UID a_ptno_Descrip_UID] = Generate_Label(apply_series(jj));
        a_title=[a_dir a_ptno_Descrip_UID '.nii.gz'];
        w_title=[w_dir a_ptno_Descrip_UID '.nii.gz'];
        
        fprintf(fid,['\t$(ANTSPATH)/antsApplyTransforms -d 3 '...
                    '-i ' a_title ' '... % input
                    '-o ' w_title ' '... % output
                    '-r ' f_title ' '... % reference to match rez
                    '-n NearestNeighbor '... % Interpolation method for resampling
                    '-t [' routput_dir 'RegistrationOutput_' m_ptno_Descrip_UID '_0GenericAffine.mat,0] '... % [..., 0 ]= Use forward transform
                    '\n']); 
    end
end

fclose(fid);
system(['make -j ' num2str(n_cores) ' -f ' prefix makefile])





end