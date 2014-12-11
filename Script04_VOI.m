%% I. Load stuff
clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';

mksqlite('open', 'ctkDICOM.sql' );
series_all = mksqlite(['select * from Series order by SeriesDate ASC']);
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);
series_SWAN = mksqlite(['select * from Series where SeriesDescription like ''%SWAN%'' order by SeriesDate ASC']);%14, none for patient 6 or 7
series_T1 = mksqlite(['select * from Series where SeriesDescription like ''%Ax T1%'' and SeriesDescription not like ''%+C%'' order by SeriesDate ASC']); %15
series_T1post = mksqlite(['select * from Series where SeriesDescription like ''%T1%'' and SeriesDescription like ''%STEALTH%'' or SeriesDescription like ''%+C%'' order by SeriesDate ASC']); %15
series_T2 = mksqlite(['select * from Series where SeriesDescription like ''%T2%'' and SeriesDescription not like ''%*%'' and SeriesDescription not like ''%FLAIR%'' order by SeriesDate ASC']); %15
series_T2star = mksqlite(['select * from Series where SeriesDescription like ''%T2*%'' order by SeriesDate ASC']); %15
series_T2FLAIR = mksqlite(['select * from Series where SeriesDescription like ''%FLAIR%'' order by SeriesDate ASC']); %15

n_series_all = size(series_all,1);
n_studies_all = size(studies_all,1);
n_SWAN = size(series_SWAN,1);
n_T1 = size(series_T1,1);
n_T1post = size(series_T1post,1);
n_T2 = size(series_T2,1);
n_T2star = size(series_T2star,1);
n_T2FLAIR = size(series_T2FLAIR,1);

% Ctrl+I - Layer Inspector with World (ITK coordinates), use mouse scroll
% wheel in all 3 views
% Paintbrush Mode 2nd-to-last button 
% Use "Center on cursor" to position points

%% II. Resample image portion at high rez
% cd /workarea/igilab/github/ExLib/vtkProjection
% git pull
% git clone http://github.com/ImageGuidedTherapyLab/ExLib.git
% chgrp igilab github/
% ll = list with details
fid = fopen([script04_prefix 'Extract.makefile'],'w');
fprintf(fid, 'all:');
for ii=1:n_series_all 
    fprintf(fid, [' job' num2str(ii)]);
end

zz=1;
for ii = 4 % patient #
    point1 = [-15.03 -53.67 58.24]; % ITK kij
    point2 = [-14.37 -52.98 46.27];
    point3 = [-13.58 -52.37 40.3]; % x y z
    start_point = point3;
    site = 'P1C';
    ptno = sprintf(['%02d'],ii);
    ptno_site = [ptno '_' site];
    
    StudyInstanceUID = studies_all(ii).StudyInstanceUID;
    temp_series = mksqlite(['select * from Series where StudyInstanceUID = ''' StudyInstanceUID ''' '...
        'and SeriesDescription like ''%FLAIR%'' '...
        'order by SeriesDate ASC']); %15
    T2FLAIR_SeriesInstanceUID = temp_series.SeriesInstanceUID; % good

    % Resample at high rez
    file = load_nii([script01_prefix '07_HistMatch/Brain_HistMatch_' T2FLAIR_SeriesInstanceUID '.nii.gz']);
    % xyz --> kij
    srow_x = file.hdr.hist.srow_x; % [ITK k]
    srow_y = file.hdr.hist.srow_y; % [ITK i]
    srow_z = file.hdr.hist.srow_z; % [ITK j]
    % ITK xyz -> ITK kij
    k = double(round(1 + (-start_point(1) - srow_x(4))/srow_x(3))); % ITK k
    i = double(round(1 + (-start_point(2) - srow_y(4))/srow_y(1))); % ITK i
    j = double(round(1 + (start_point(3) - srow_z(4))/srow_z(2))); % ITK j
    extract_side = 12; %mm
    origin_k = k-(extract_side/abs(srow_x(3)))/2;
    origin_i = i-(extract_side/abs(srow_y(1)))/2;
    origin_j = j-(extract_side/abs(srow_z(2)))/2;
    num_k = extract_side/abs(srow_x(3));
    num_i = extract_side/abs(srow_y(1));
    num_j = extract_side/abs(srow_z(2));
    extract_origin = [num2str(origin_i) 'x' num2str(origin_j) 'x' num2str(origin_k)];
    extract_size = [num2str(num_i) 'x' num2str(num_j) 'x' num2str(num_k)];
    
    % T2FLAIR
    T2FLAIR_SeriesDescription = strrep(strrep(strrep(strrep(strrep(temp_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    T2FLAIR_file_infix = sprintf([ T2FLAIR_SeriesDescription '_' T2FLAIR_SeriesInstanceUID]); % patient #, Descrip, SeriesInstanceUID
    fprintf(fid,['\n\njob' num2str(zz) ':\n']);
    fprintf(fid,['\tc3d ' script01_prefix '07_HistMatch/Brain_HistMatch_' T2FLAIR_SeriesInstanceUID '.nii.gz '... % Input filename
    		'-region ' extract_origin 'vox ' extract_size 'vox '... % Origin-ijk dimensions-ijk
            '-type double -noround '...
    		'-o ' script01_prefix '11_Extract/Extract_' ptno_site '_' T2FLAIR_file_infix '.nii.gz\n']);
    fprintf(fid,['\tc3d ' script01_prefix '11_Extract/Extract_' ptno_site '_' T2FLAIR_file_infix '.nii.gz '... % Input filename
    		'-interpolation NearestNeighbor -resample-mm 0.02x0.02x0.02mm '... % Interpolation type, resampling dimensions
    		'-type double -noround '...
            '-o ' script01_prefix '12_Resample/Resample_' ptno_site '_' T2FLAIR_file_infix '.nii.gz\n']);
        
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script01_prefix '10_Warped/Warped_' file_infix '.nii.gz '])
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script01_prefix '11_Extract/Extract_' ptno_site '_' file_infix '.nii.gz'])
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script01_prefix '12_Resample/Resample_' ptno_site '_' T2FLAIR_file_infix '.nii.gz']) 
        
        
        zz=zz+1;

    % T1,T2,SWAN
    temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%SWAN%'' '... % Includes SWAN Image volume
        'or SeriesDescription like ''%T1%'' '... % Includes T1, T1post
        'or SeriesDescription like ''%T2%'') '... % Includes T2, T2star, T2FLAIR
        'and SeriesDescription not like ''%FLAIR%'' '... % Exclude FLAIR - different save location
        'order by SeriesNumber ASC']);
    for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
        UID = temp_series(jj).SeriesInstanceUID;
        Descrip = strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
        Descrip_UID = sprintf([Descrip '_' UID]); % Descrip, SeriesInstanceUID
        
        fprintf(fid,['\n\njob' num2str(zz) ':\n']);
        fprintf(fid,['\tc3d ' script01_prefix '10_Warped/Warped_' Descrip_UID '.nii.gz '... % Input filename
                '-region ' extract_origin 'vox ' extract_size 'vox '... % Origin-ijk dimensions-ijk
                '-type double -noround '...
                '-o ' script01_prefix '11_Extract/Extract_' ptno_site '_' Descrip_UID '.nii.gz\n']);
        fprintf(fid,['\tc3d ' script01_prefix '11_Extract/Extract_' ptno_site '_' Descrip_UID '.nii.gz '... % Input filename
                '-interpolation NearestNeighbor -resample-mm 0.02x0.02x0.02mm '... % Interpolation type, resampling dimensions
                '-type double -noround '...
                '-o ' script01_prefix '12_Resample/Resample_' ptno_site '_' Descrip_UID '.nii.gz\n']);
            
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script01_prefix '07_HistMatch/Brain_HistMatch_' T2FLAIR_SeriesInstanceUID '.nii.gz '...
%             '-o ' script01_prefix '10_Warped/Warped_' file_infix '.nii.gz ' ]) % Check Registration
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script01_prefix '10_Warped/Warped_' file_infix '.nii.gz '])
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script01_prefix '11_Extract/Extract_' ptno_site '_' file_infix '.nii.gz'])
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script01_prefix '12_Resample/Resample_' ptno_site '_' file_infix '.nii.gz']) 
        
        zz=zz+1;
    end

    % DWI/DTI
    temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
            'and (SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
            'or SeriesDescription like ''%Fractional Aniso%'' '... % Includes FA from DTI
            'or SeriesDescription like ''%Average DC%'') '...  % Includes ADC from DTI
            'order by SeriesNumber ASC']);
    for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
        UID = temp_series(jj).SeriesInstanceUID;
        Descrip = strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    
        Descrip_UID = sprintf([Descrip '_' UID]);
        fprintf(fid,['\n\njob' num2str(zz) ':\n']);
        fprintf(fid,['\tc3d ' script02_prefix '10_Warped/Warped_' ptno '_' Descrip_UID '.nii.gz '... % Input filename
                '-region ' extract_origin 'vox ' extract_size 'vox '... % Origin-ijk dimensions-ijk
                '-type double -noround '...
                '-o ' script02_prefix '11_Extract/Extract_' ptno_site '_' Descrip_UID '.nii.gz\n']);
        fprintf(fid,['\tc3d ' script02_prefix '11_Extract/Extract_' ptno_site '_' Descrip_UID '.nii.gz '... % Input filename
                '-interpolation NearestNeighbor -resample-mm 0.02x0.02x0.02mm '... % Interpolation type, resampling dimensions
                '-type double -noround '...
                '-o ' script02_prefix '12_Resample/Resample_' ptno_site '_' Descrip_UID '.nii.gz\n']);
            
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script01_prefix '07_HistMatch/Brain_HistMatch_' T2FLAIR_SeriesInstanceUID '.nii.gz '...
%             '-o ' script02_prefix '10_Warped/Warped_' ptno '_' file_infix '.nii.gz ' ]) % Check Registration
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script02_prefix '10_Warped/Warped_' ptno '_' file_infix '.nii.gz'])
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script02_prefix '11_Extract/Extract_' ptno_site '_' file_infix '.nii.gz'])
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script02_prefix '12_Resample/Resample_' ptno_site '_' file_infix '.nii.gz']) 
        
        zz=zz+1;
    end

    
    % DSC/DCE
    temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%rBV%'' '... % DSC maps
        'or SeriesDescription like ''%rBF%'' '...
        'or SeriesDescription like ''%MTT%'' '...
        'or SeriesDescription like ''%Delay%'' '...
        'or SeriesDescription like ''%Leakage (K2)%'' '...
        'or SeriesDescription like ''%K12%'' '... % DCE maps
    	'or SeriesDescription like ''%K21%'' '...
        'or SeriesDescription like ''%Vp%'' '...
    	'or SeriesDescription like ''%Ve (distribution volume)%'' '...
    	'or SeriesDescription like ''%Wash-in%'' '...
    	'or SeriesDescription like ''%Wash-out%'' '...
    	'or SeriesDescription like ''%Time to peak%'' '...
    	'or SeriesDescription like ''%Area under curve%'' '...
    	'or SeriesDescription like ''%Peak enhancement%'')'...
        'order by SeriesNumber ASC']);
    for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
        UID = temp_series(jj).SeriesInstanceUID;
        Descrip = strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
    
        Descrip_UID = sprintf([Descrip '_' UID]);
        fprintf(fid,['\n\njob' num2str(zz) ':\n']);
        fprintf(fid,['\tc3d ' script03_prefix '10_Warped/Warped_' ptno '_' Descrip_UID '.nii.gz '... % Input filename
                '-region ' extract_origin 'vox ' extract_size 'vox '... % Origin-ijk dimensions-ijk
                '-type double -noround '...
                '-o ' script03_prefix '11_Extract/Extract_' ptno_site '_' Descrip_UID '.nii.gz\n']);
        fprintf(fid,['\tc3d ' script03_prefix '11_Extract/Extract_' ptno_site '_' Descrip_UID '.nii.gz '... % Input filename
                '-interpolation NearestNeighbor -resample-mm 0.02x0.02x0.02mm '... % Interpolation type, resampling dimensions
                '-type double -noround '...
                '-o ' script03_prefix '12_Resample/Resample_' ptno_site '_' Descrip_UID '.nii.gz\n']);
            
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script01_prefix '07_HistMatch/Brain_HistMatch_' T2FLAIR_SeriesInstanceUID '.nii.gz '...
%             '-o ' script03_prefix '10_Warped/Warped_' ptno '_' file_infix '.nii.gz ']) % Check Registration
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script03_prefix '10_Warped/Warped_' ptno '_' file_infix '.nii.gz'])
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script03_prefix '11_Extract/Extract_' ptno_site '_' file_infix '.nii.gz'])
%         system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%             '-g ' script03_prefix '12_Resample/Resample_' ptno_site '_' file_infix '.nii.gz'])     
            
        zz=zz+1;
    end


end

fclose(fid);
system(['make -j 8 -f ' script04_prefix 'Extract.makefile'])
% 3.88 min

system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '07_HistMatch/Brain_HistMatch_' T2FLAIR_SeriesInstanceUID '.nii.gz '])
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script04_prefix 'Brain_HistMatch_' T2FLAIR_SeriesInstanceUID '_Extract.nii.gz'])
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script04_prefix 'Brain_HistMatch_' T2FLAIR_SeriesInstanceUID '_Resample.nii.gz']) % 0.6 min













%% III. Create VOI

ii = 4; % patient #
point1 = [-15.03 -53.67 58.24]; % ITK kij
point2 = [-14.37 -52.98 46.27];
point3 = [-13.58 -52.37 40.3]; % x y z
start_point = point3;

% Get trajectory line based off 2 points
ctr_vector1 = point1-point2;
ctr_vector2 = point2-point3;
ctr_vector3 = point1-point3;
ctr_vector_unit = ctr_vector2/norm(ctr_vector2);

% Create labels
site = 'P1C';
ptno = sprintf(['%02d'],ii);
ptno_site = sprintf(['%02d_' site],ii);
ptno_site_mask = [ptno_site '_Mask'];

temp_series = mksqlite(['select * from Series where StudyInstanceUID = ''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%FLAIR%'' '...
    'order by SeriesDate ASC']); %15
UID = temp_series.SeriesInstanceUID; % good
Descrip = strrep(strrep(strrep(strrep(strrep(temp_series.SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star');
Descrip_UID = sprintf([ Descrip '_' UID]); % patient #, Descrip, SeriesInstanceUID


for rez = 1 % rez 1 Original, 2 Resampled
% Resample at high rez

if rez==1
    file = load_nii([script01_prefix '07_HistMatch/Brain_HistMatch_' UID '.nii.gz']);
elseif rez==2
    file = load_nii([script01_prefix '12_Resample/Resample_' ptno_site '_' Descrip_UID '.nii.gz']);
end

imgdim = [file.hdr.dime.dim(2:4)]; % # pixels - ijk
srow_x = file.hdr.hist.srow_x; % [ITK k] % xyz --> kij
srow_y = file.hdr.hist.srow_y; % [ITK i]
srow_z = file.hdr.hist.srow_z; % [ITK j]
mask = zeros(imgdim);

choice = 1; % cylinder = 1, Forceps = 2
if choice == 1 % Cylinder
    mask_cyl = mask;
    tt_ctr_limit = 4.5; % height/2 = 9/2 = 4.5 mm
    step = 0.01; % mm
    for tt_ctr  = -tt_ctr_limit:step:tt_ctr_limit % inch along cylinder
        tt_ctr
        ctr_point = start_point + ctr_vector_unit*tt_ctr;
        a = ctr_point(1); % Point on line
        b = ctr_point(2);
        c = ctr_point(3);

        % Central axis vector
        u = ctr_vector_unit(1); 
        v = ctr_vector_unit(2);
        w = ctr_vector_unit(3);

	% Perpendicular vector = [1 1 -(A+B)/C]
        x = 1; 
        y = 1;
        z = -(u + v)/w;

        deg = 0:1:360; % degrees by which it should be rotated
        ct = cos(deg*pi/180);
        st = sin(deg*pi/180);

        % Rotate 360 degrees, simultaneously
        edge_vector_rot = [ (a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-ct)+x*ct+(-c*v+b*w-w*y+v*z)*st
        (b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-ct)+y*ct+(c*u-a*w+w*x-u*z)*st
        (c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-ct)+z*ct+(-b*u+a*v-v*x+u*y)*st ]';
        norm_evr=sqrt(sum(edge_vector_rot.^2,2));
        evr_unit = edge_vector_rot./[norm_evr norm_evr norm_evr];
        edge_limit = .735; % radius = 1.47/2 = .7350

        for tt_edge = -edge_limit:step:edge_limit 
            edge_point = repmat([a b c],size(evr_unit,1),1) + evr_unit*tt_edge;
            x_new = edge_point(:,1); 
            y_new = edge_point(:,2); 
            z_new = edge_point(:,3); 

            % ITK xyz -> ITK kij
            k = double(round(1 + (-x_new - srow_x(4))/srow_x(3))); % ITK k
            i = double(round(1 + (-y_new - srow_y(4))/srow_y(1))); % ITK i
            j = double(round(1 + (z_new - srow_z(4))/srow_z(2))); % ITK j
            idx_cyl = sub2ind(size(mask_cyl),i,j,k);
            mask_cyl(idx_cyl)=1; % for close voxels, mask = 1

        end % inching along edge vector
    end % inch along center line

    holder_cyl = make_nii(mask,... % Image
        [file.hdr.dime.pixdim(2:4)], ... % Pixel dimensions - ijk
        [srow_x(4) srow_y(4) srow_z(4)],... % origin - xyz
        64,[ptno_site_mask '_Cylinder']); % datatype: 64 = float64 = double, Description
    holder_cyl.hdr = file.hdr;

elseif choice == 2 % Forceps
    mask_forceps = mask;
    [j_all i_all k_all] = meshgrid(1:imgdim(2),1:imgdim(1),1:imgdim(3));
    
    % IJK coordinates --> XYZ
    x_all = -((srow_x(1)*(i_all-1)) + (srow_x(2)*(j_all-1)) + (srow_x(3)*(k_all-1)) + srow_x(4));
    y_all = -((srow_y(1)*(i_all-1)) + (srow_y(2)*(j_all-1)) + (srow_y(3)*(k_all-1)) + srow_y(4));
    z_all = (srow_z(1)*(i_all-1)) + (srow_z(2)*(j_all-1)) + (srow_z(3)*(k_all-1)) + srow_z(4);
    
    % Calculate XYZ coordinates, check if < 5mm distance from start_point
    x_diff = x_all - start_point(1);
    y_diff = y_all - start_point(2);
    z_diff = z_all - start_point(3);
    xyz_dist = sqrt(x_diff.^2+y_diff.^2+z_diff.^2);   
    idx_forceps = find(xyz_dist(:) < .5); % .5 mm radius, 1 mm diameter
    mask_forceps(idx_sphere)=1; % for close voxels, mask = 1

    holder_forceps = make_nii(mask,... % Image
        [file.hdr.dime.pixdim(2:4)], ... % Pixel dimensions - ijk
        [srow_x(4) srow_y(4) srow_z(4)],... % origin - xyz
        64,[ptno_site_mask '_Forceps']); % datatype: 64 = float64 = double, Description
    holder_forceps.hdr = file.hdr;
    
end

    % Sphere - performed in every case
    mask_sphere = mask;
    [j_all i_all k_all] = meshgrid(1:imgdim(2),1:imgdim(1),1:imgdim(3));
    
    % IJK coordinates --> XYZ
    x_all = -((srow_x(1)*(i_all-1)) + (srow_x(2)*(j_all-1)) + (srow_x(3)*(k_all-1)) + srow_x(4));
    y_all = -((srow_y(1)*(i_all-1)) + (srow_y(2)*(j_all-1)) + (srow_y(3)*(k_all-1)) + srow_y(4));
    z_all = (srow_z(1)*(i_all-1)) + (srow_z(2)*(j_all-1)) + (srow_z(3)*(k_all-1)) + srow_z(4);
    
    % Calculate XYZ coordinates, check if < 5mm distance from start_point
    x_diff = x_all - start_point(1);
    y_diff = y_all - start_point(2);
    z_diff = z_all - start_point(3);
    xyz_dist = sqrt(x_diff.^2+y_diff.^2+z_diff.^2);   
    idx_sphere = find(xyz_dist(:) < 2.5); % 2.5 mm radius, 5 mm diameter
    mask_sphere(idx_sphere)=1; % for close voxels, mask = 1
    
    holder_sphere = make_nii(mask,... % Image
        [file.hdr.dime.pixdim(2:4)], ... % Pixel dimensions - ijk
        [srow_x(4) srow_y(4) srow_z(4)],... % origin - xyz
        64,[ptno_site_mask '_Sphere']); % datatype: 64 = float64 = double, Description
    holder_sphere.hdr = file.hdr; % Valid bc I commented out the x_form_nii command in load_nii,which changed qform/sform and applied a transform to the data

    % Commented out the hdr qform, sform assignments in save_nii_hdr so it would stop being annoying
    script04_prefix_results = [script04_prefix 'Results/' ptno '/' ptno_site '/'];
    if ~isdir(script04_prefix_results)
        mkdir(script04_prefix_results)
    end

    if choice==1 && rez==1 % Cylinder/Sphere, Original Rez
        save_nii(holder_cyl, [script04_prefix_results ptno_site_mask '_Cyl_Original.nii.gz']) % 2.15 min - Commented out the hdr qform, sform
        save_nii(holder_sphere, [script04_prefix_results ptno_site_mask '_Sphere_Original.nii.gz']) % 2.15 min - Commented out the hdr qform, sform
    elseif choice==1 && rez==2 % Cylinder/Sphere, Resampled Rez
        save_nii(holder_cyl, [script04_prefix_results ptno_site_mask '_Cyl_Resample.nii.gz']) 
        save_nii(holder_sphere, [script04_prefix_results ptno_site_mask '_Sphere_Resample.nii.gz']) 
    elseif choice==2 && rez==1 % Forceps/Sphere, Original Rez
        save_nii(holder_forceps, [script04_prefix_results ptno_site_mask '_Forceps_Original.nii.gz']) 
        save_nii(holder_sphere, [script04_prefix_results ptno_site_mask '_Sphere_Original.nii.gz']) 
    elseif choice==2 && rez==2 % Forceps/Sphere, Resampled Rez
        save_nii(holder_forceps, [script04_prefix_results ptno_site_mask '_Forceps_Resample.nii.gz']) 
        save_nii(holder_sphere, [script04_prefix_results ptno_site_mask '_Sphere_Resample.nii.gz']) 
    end


end % Perform at rez = 1 (Original), then rez = 2 (Resampled)

    %%
    system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
        '-g ' script04_prefix ptno_site_mask '_Cyl.nii.gz']) % 1.15 min
    system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
        '-g ' file.fileprefix '.nii.gz '...
        '-s ' script04_prefix ptno_site_mask '_Cyl.nii.gz']) % 1.15 min
    
    system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
        '-g ' script04_prefix_results ptno_site_mask '_Sphere.nii.gz']) % 1.15 min
    system(['vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap '...
        '-g ' file.fileprefix '.nii.gz '...
        '-s ' script04_prefix_results ptno_site_mask '_Sphere.nii.gz']) % 1.15 min

close all
% Verified same coordinates in Matlab and ITKsnap format - 20141105

for qq = 291%[260 280 305 313 340] %coronal - 315:5:365
    figure; imagesc(squeeze(file.img(qq,:,:))); colormap gray; axis off; title(num2str(qq))
    figure; imagesc(squeeze(mask(qq,:,:))); colormap gray; axis off; title(num2str(qq))
%     figure; imagesc(squeeze(img_x_mask_resample.img(qq,:,:))); colormap gray; axis off; title(num2str(qq))
end

for ww = 116%[313 340 366 392]% axial - 315:5:335
    figure; imagesc(squeeze(file.img(:,ww,:))); colormap gray; axis off; title(num2str(ww))
    figure; imagesc(squeeze(mask(:,ww,:))); colormap gray; axis off; title(num2str(ww))
%     figure; imagesc(squeeze(img_x_mask_resample.img(:,ww,:))); colormap gray; axis off; title(num2str(ww))
end

for ee = 93%[290 300 313 330 350]% sagittal - 265:5:290
    figure; imagesc(squeeze(file.img(:,:,ee))); colormap gray; axis off; title(num2str(ee))
    figure; imagesc(squeeze(mask(:,:,ee))); colormap gray; axis off; title(num2str(ee))
%     figure; imagesc(squeeze(img_x_mask_resample.img(:,:,ee))); colormap gray; axis off; title(num2str(ee))
end






%% V. Use VOI mask to save mean, sd, histogram in .mat file


tic
for rez = 1:2

    zz=1;
    voi_cell1 = cell(1,9);
    voi_cell2 = cell(1,9);

    if rez == 1
        if choice == 1
            mask = load_nii([script04_prefix_results ptno_site_mask '_Cyl_Original.nii.gz']); 
        elseif choice == 2
            mask = load_nii([script04_prefix_results ptno_site_mask '_Forceps_Original.nii.gz']); 
        end
        mask_sphere = load_nii([script04_prefix_results ptno_site_mask '_Sphere_Original.nii.gz']); 
    elseif rez == 2
        if choice == 1
            mask = load_nii([script04_prefix_results ptno_site_mask '_Cyl_Resample.nii.gz']); 
        elseif choice == 2
            mask = load_nii([script04_prefix_results ptno_site_mask '_Forceps_Resample.nii.gz']); 
        end
        mask_sphere = load_nii([script04_prefix_results ptno_site_mask '_Sphere_Resample.nii.gz']); 
    end

voi_idx1=find(mask.img);
voi_idx2=find(mask_sphere.img);



% T1
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%T1%'' '... % Includes T1, T1post
    'and (SeriesDescription not like ''%STEALTH%'' '... %Excludes post
    'and SeriesDescription not like ''%+C%'') '... % Excludes post
    'order by SeriesNumber ASC']);
VOI_output1(temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
zz = zz+1;

% T1post
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%T1%'' '... % Includes T1, T1post
    'and (SeriesDescription like ''%STEALTH%'' '... %Includes only post
    'or SeriesDescription like ''%+C%'') '... % Includes only post
    'order by SeriesNumber ASC']);
VOI_output1(temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
zz = zz+1;

% T2
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%T2%'' '... % Includes SWAN Image volume
    'and SeriesDescription not like ''%*%'' '... % Includes T1, T1post
    'and SeriesDescription not like ''%FLAIR%'' '... % Includes T2, T2star, T2FLAIR
    'order by SeriesNumber ASC']);
VOI_output1(temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
zz = zz+1;

% T2star
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%T2*%'' '... % Includes SWAN Image volume
    'order by SeriesNumber ASC']);
VOI_output1(temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
zz = zz+1;

% T2 FLAIR
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%FLAIR%'' '... % Includes T2, T2star, T2FLAIR
    'order by SeriesNumber ASC']);
VOI_output1(temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
zz = zz+1;

% SWAN
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%SWAN%'' '... % Includes SWAN Image volume
    'order by SeriesNumber ASC']);
VOI_output1(temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
zz = zz+1;



%% DWI/DTI
conv_factor = 1e-6;

%ADC
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
        'and SeriesDescription not like ''%Exponential%'' '... % Includes FA from DTI
        'order by SeriesNumber ASC']);
VOI_output2(1e-6, temp_series, rez, ptno, site,  voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
zz = zz+1;
   
% eADC    
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%Exponential Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
        'order by SeriesNumber ASC']);   
VOI_output2(1e-6, temp_series, rez, ptno, site,  voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
zz = zz+1;
 
% FA
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%Fractional Aniso%'' '... % Includes FA from DTI
        'order by SeriesNumber ASC']);   
VOI_output2(1e-6, temp_series, rez, ptno, site,  voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
zz = zz+1;
 
% Avg DC    
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and SeriesDescription like ''%Average DC%'' '...  % Includes ADC from DTI
        'order by SeriesNumber ASC']);   
VOI_output2(1e-6, temp_series, rez, ptno, site,  voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
zz = zz+1;

    
%% DSC
%rBV - not corrected
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%rBV%'' '... % DSC maps
    'and SeriesDescription like ''%not%'' '... % DSC maps
    'order by SeriesNumber ASC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

%rBV_corr
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%rBV%'' '... % DSC maps
    'and SeriesDescription not like ''%not%'' '... % DSC maps
    'order by SeriesNumber ASC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

%rBF
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%rBF%'' '...
    'order by SeriesNumber ASC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% MTT
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%MTT%'' '...
    'order by SeriesNumber ASC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% Delay
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%Delay%'' '...
    'order by SeriesNumber ASC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% K2
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%Leakage (K2)%'' '...
    'order by SeriesNumber ASC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

%% DCE
% Ktrans
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%K12%'' '... % DCE maps
    'order by SeriesDescription DESC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% Kep
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%K21%'' '...
    'order by SeriesDescription DESC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% Vp
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%Vp%'' '...
    'order by SeriesDescription DESC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% Ve
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%Ve (distribution volume)%'' '...
    'order by SeriesDescription DESC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% Wash-in
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%Wash-in%'' '...
    'order by SeriesDescription DESC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% Wash-out
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%Wash-out%'' '...
    'order by SeriesDescription DESC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% TTP
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%Time to peak%'' '...
    'order by SeriesDescription DESC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% AUC
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%Area under curve%'' '...
    'order by SeriesDescription DESC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

% Peak enhancement
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
    'and SeriesDescription like ''%Peak enhancement%'' '...
    'order by SeriesDescription DESC']);
for jj = 1:size(temp_series,1) % will atuomatically pass over if no content in temp_series
    zz
    VOI_output3(jj, temp_series, rez, ptno, site, voi_cell1, voi_cell2, voi_idx1, voi_idx2, zz)
    zz = zz+1;
end

    % Save voi_cell1 and voi_cell2 as .mat files
    if choice==1 && rez==1 % Cylinder/Sphere, Original Rez
        save([script04_prefix_results ptno '_' site '_Cyl_Original.mat'], 'voi_cell1') 
        save([script04_prefix_results ptno '_' site '_Sphere_Original.mat'], 'voi_cell2') 
    elseif choice==1 && rez==2 % Cylinder/Sphere, Resampled Rez
        save([script04_prefix_results ptno '_' site '_Cyl_Resample.mat'], 'voi_cell1') 
        save([script04_prefix_results ptno '_' site '_Sphere_Resample.mat'], 'voi_cell2') 
    elseif choice==2 && rez==1 % Forceps/Sphere, Original Rez
        save([script04_prefix_results ptno '_' site '_Forceps_Original.mat'], 'voi_cell1') 
        save([script04_prefix_results ptno '_' site '_Sphere_Original.mat'], 'voi_cell2') 
    elseif choice==2 && rez==2 % Forceps/Sphere, Resampled Rez
        save([script04_prefix_results ptno '_' site '_Forceps_Resample.mat'], 'voi_cell1') 
        save([script04_prefix_results ptno '_' site '_Sphere_Resample.mat'], 'voi_cell2') 
    end

end % end rez outer loop


toc/60 % 26.74 min

close all

% sort into this order
% 1. T1,
% 2. T1post,
% 3. T2,
% 4. T2star,
% 5. T2FLAIR,
% 6. SWAN
% 7. ADC, % DWI
% 8. eADC, 
% 9. FA, % DTI
% 10. AvgDC
% 11. rBV, % DSC
% 12. rBF, 
% 13. MTT, 
% 14. Delay, 
% 15. Leakage, 
% 16. K12,  % DCE
% 17. K21, 
% 18. Vp, 
% 19. Ve, 
% 20. Wash-in, 
% 21. Wash-out, 
% 22. TTP, 
% 23. AUC, 
% 24. Peak Enhancement

% Check if registrations went okay?


%% IV. Check VOI on top of each image volume
% Slice through mask to see if orientation is same as resampled extract

clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/mksqlite-1.14')
path(path, '/mnt/data/scratch/igilab/jslin1/Matlab_add-ons/NIfTI_20140122')
script01_prefix = 'Script01_T1_T2_SWAN/';
script02_prefix = 'Script02_DWI_DTI/';
script03_prefix = 'Script03_DCE_DSC/';
script04_prefix = 'Script04_VOI/';

mksqlite('open', 'ctkDICOM.sql' );
series_all = mksqlite(['select * from Series order by SeriesDate ASC']);
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);

ii=4;
site = 'P1C';
ptno = sprintf(['%02d'],ii);
ptno_site = sprintf(['%02d_' site],ii);
ptno_site_mask = [ptno_site '_Mask'];

% Anatomical
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%T1%'' '... % Includes T1, T1post
        'or SeriesDescription like ''%T2%'' '... % Includes T2, T2star, T2FLAIR
        'or SeriesDescription like ''%SWAN%'')'...
        'order by SeriesDate, SeriesInstanceUID ASC']);
jj=4;
UID = temp_series(jj).SeriesInstanceUID; % good
Descrip = strrep(strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star'),'.','');
Descrip_UID = sprintf([ Descrip '_' UID]); % patient #, Descrip, SeriesInstanceUID
% file2 = load_nii([script01_prefix '10_Warped/Warped_' Descrip_UID '.nii.gz']);
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '10_Warped/Warped_' Descrip_UID '.nii.gz '...
    '-s ' script04_prefix ptno_site_mask '.nii.gz']) % 1.15 min
% system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
%     '-g ' script01_prefix '07_HistMatch/Brain_HistMatch_' UID '.nii.gz '...
%     '-s ' script04_prefix mask_title '.nii.gz']) % 1.15 min

system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script01_prefix '12_Resample/Resample_' ptno_site '_' Descrip_UID '.nii.gz '...
    '-s ' script04_prefix ptno_site_mask '_Resample.nii.gz']) % 1.15 min
% Bad registration - T1
% 4.8X

% DWI/DTI
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%Apparent Diffusion Coefficient%'' '... % Includes ADC and eADC from DWI
        'or SeriesDescription like ''%Fractional Aniso%'' '... % Includes FA from DTI
        'or SeriesDescription like ''%Average DC%'')'...
        'order by SeriesDate, SeriesInstanceUID ASC']);
jj=4;
UID = temp_series(jj).SeriesInstanceUID; % good
Descrip = strrep(strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star'),'.','');
Descrip_UID = sprintf([ Descrip '_' UID]); % patient #, Descrip, SeriesInstanceUID
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script02_prefix '10_Warped/Warped_' ptno '_' Descrip_UID '.nii.gz '...
    '-s ' script04_prefix ptno_site_mask '.nii.gz']) % 1.15 min
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script02_prefix '12_Resample/Resample_' ptno_site '_' Descrip_UID '.nii.gz '...
    '-s ' script04_prefix ptno_site_mask '_Resample.nii.gz']) % 1.15 min

% DSC/DCE
temp_series = mksqlite(['select * from Series where StudyInstanceUID=''' studies_all(ii).StudyInstanceUID ''' '...
        'and (SeriesDescription like ''%DCE%'' '...
        'or SeriesDescription like ''%DSC%'' '...
        'or SeriesDescription like ''%rBV%'' '...
        'or SeriesDescription like ''%rBF%'' '...
        'or SeriesDescription like ''%MTT%'' '...
        'or SeriesDescription like ''%Delay%'' '...
        'or SeriesDescription like ''%Leakage (K2)%'' '...
        'or SeriesDescription like ''%K12%'' '...
        'or SeriesDescription like ''%K21%'' '...
        'or SeriesDescription like ''%Vp%'' '...
        'or SeriesDescription like ''%Ve (distribution volume)%'' '...
        'or SeriesDescription like ''%Wash-in%'' '...
        'or SeriesDescription like ''%Wash-out%'' '...
        'or SeriesDescription like ''%Time to peak%'' '...
        'or SeriesDescription like ''%Area under curve%'' '...
        'or SeriesDescription like ''%Peak enhancement%'')'...
        'order by SeriesDate, SeriesInstanceUID ASC']);
jj=3;
UID = temp_series(jj).SeriesInstanceUID; % good
Descrip = strrep(strrep(strrep(strrep(strrep(strrep(temp_series(jj).SeriesDescription,' ','_'),'(','_'),')','_'),'/','_'),'*','star'),'.','');
Descrip_UID = sprintf([ Descrip '_' UID]); % patient #, Descrip, SeriesInstanceUID
system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script03_prefix '10_Warped/Warped_' ptno '_' Descrip_UID '.nii.gz '...
    '-s ' script04_prefix ptno_site_mask '.nii.gz']) % 1.15 min

system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script03_prefix '10_Warped/Warped_' ptno '_' Descrip '_lasttimept_' UID '.nii.gz '...
    '-s ' script04_prefix ptno_site_mask '.nii.gz']) % 1.15 min

system(['vglrun /opt/apps/itksnap/itksnap-3.0.0-20140425-Linux-x86_64/bin/itksnap '...
    '-g ' script03_prefix '12_Resample/Resample_' ptno_site '_' Descrip_UID '.nii.gz '...
    '-s ' script04_prefix ptno_site_mask '_Resample.nii.gz']) % 1.15 min





