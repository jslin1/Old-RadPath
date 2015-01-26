function mask_cyl = Mask_Cyl(ii, point_info, series_sql, x_all, x_step,y_step,z_step,side)

path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql

%% ptno --> Correct sql entry
studies_all = mksqlite(['select * from Studies order by StudyDate ASC']);

% Ptno -> StudyInstanceUID -> series_sql entry
ptno = cell2mat(point_info(ii,1));
% series_sql = series_T2_Register; % for troubleshooting
[field study]= find(cell2mat(cellfun(@(x) isequal(studies_all(str2num(ptno)).StudyInstanceUID,x),struct2cell(series_sql),'UniformOutput',false)));
[ptno Descrip UID ptno_Descrip_UID] = Generate_Label(series_sql(study));

% Load srow_xyz
file = load_nii([ script01_prefix '01_N4/N4_' ptno_Descrip_UID '.nii.gz' ]);
srow_x = file.hdr.hist.srow_x; % [ITK k]
srow_y = file.hdr.hist.srow_y; % [ITK i]
srow_z = file.hdr.hist.srow_z; % [ITK j]

% Point and vcp info
point = cell2mat(point_info(ii,4:6));
vcp = cell2mat(point_info(ii,7));

%% Define Center Vector using vcp (vector control point)
if isnumeric(vcp) % Tip Coordinates
    vcp_num = cell2mat(point_info(ii,7:9));
elseif ischar(vcp) % '01-P1C'
    vcp_char = vcp;
    idx_row=find(cell2mat(cellfun(@(x) isequal(vcp_char,x),point_info(:,2),'UniformOutput',false)));
    vcp_num = cell2mat(point_info(idx_row,4:6)); % Find/oad correct row
end
ctr_vector = point-vcp_num;
ctr_vector_unit = ctr_vector/sqrt(sum(ctr_vector.^2,2));

%% Needle/Cylinder
mask_cyl = zeros(size(x_all)); % .02 Vox mask
tt_ctr_limit = 4.5; % height/2 = 9mm/2 = 4.5 and NOT 8.38 mm/2 = 4.19 mm
cyl_step = 0.01; % mm
for tt_ctr  = -tt_ctr_limit:cyl_step:tt_ctr_limit % inch along cylinder
    tt_ctr;
    ctr_point = point + ctr_vector_unit*tt_ctr;
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

    deg = 0:.8:360; % degrees by which it should be rotated
    ct = cos(deg*pi/180);
    st = sin(deg*pi/180);

    % Rotate 360 degrees, simultaneously
    edge_vector_rot = [ (a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-ct)+x*ct+(-c*v+b*w-w*y+v*z)*st
    (b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-ct)+y*ct+(c*u-a*w+w*x-u*z)*st
    (c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-ct)+z*ct+(-b*u+a*v-v*x+u*y)*st ]';
    norm_evr=sqrt(sum(edge_vector_rot.^2,2));
    evr_unit = edge_vector_rot./[norm_evr norm_evr norm_evr];
    edge_limit = .735; % radius = 1.47mm/2 = .735 mm

    for tt_edge = -edge_limit:cyl_step:edge_limit 
        edge_point = repmat([a b c],size(evr_unit,1),1) + evr_unit*tt_edge;
        x_new = edge_point(:,1); 
        y_new = edge_point(:,2); 
        z_new = edge_point(:,3); 

        %% Random World -> .02 Vox 
        % Construct own key for .02 Vox meshgrid  -
        % Why does this work? Do I have to find out? Maybe later
        srow_x_temp = [sign(srow_x(1:3))*x_step (point(1)-(side/2))*sign(srow_x(srow_x(1:3)~=0))];
        srow_y_temp = [sign(srow_y(1:3))*y_step (point(2)-(side/2))*sign(srow_y(srow_y(1:3)~=0))];
        srow_z_temp = [sign(srow_z(1:3))*z_step (point(3)-(side/2))*sign(srow_z(srow_z(1:3)~=0))];

        [i_temp j_temp k_temp] = World_to_Vox(x_new, y_new, z_new, srow_x_temp, srow_y_temp, srow_z_temp); 

        idx_temp = sub2ind(size(mask_cyl),i_temp,j_temp,k_temp);
        mask_cyl(idx_temp)=1; % mask = 1 for cylinder vox

    end % inching along edge vector
end % inch along center line
    
end