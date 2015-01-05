function mask_cyl = Mask_Cyl(ii, point_info, x_all)

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
ctr_vector_unit = ctr_vector/norm(ctr_vector);

%% Needle/Cylinder
mask_cyl = zeros(size(x_all)); % .02 Vox mask
tt_ctr_limit = 0%4.19; % height/2 = 8.38 mm/2 = 4.19 mm
step = 0.01; % mm
for tt_ctr  = -tt_ctr_limit:step:tt_ctr_limit % inch along cylinder
    tt_ctr
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
    edge_limit = .1%.675; % radius = 1.35 mm/2 = 0.675 mm

    for tt_edge = -edge_limit:step:edge_limit 
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