dim = 500;
mask = zeros(dim,150,dim);

vector = [1 0 1];
test_unit = vector/norm(vector);
u = test_unit(1);%vector_line_unit(1); % unit vector of line
v = test_unit(2);%vector_line_unit(2);
w = test_unit(3);%vector_line_unit(3);
limit1 = 100;

for qq = 0:1:limit1 % inch along cylinder
    ctrpoint = [round(dim/2) 75 round(dim/2)] + test_unit*qq;
    a = ctrpoint(1); % Point on line
    b = ctrpoint(2);
    c = ctrpoint(3);

    % close all
    % figure
    % quiver3(a,b,c,u,v,w)
    % xlim([0 dim])
    % ylim([0 dim])
    % zlim([0 dim])
    % hold
    % quiver3(a,b,c,x,y,z)

    for dd = 1:1:360 % rotate vector
        vector_perp1 = [1 1 -(vector(2) + vector(1))/vector(3)];
        x = vector_perp1(1); % vector to circumference point
        y = vector_perp1(2);
        z = vector_perp1(3);
        deg = dd; % degrees by which it should be rotated
        ct = cos(deg*pi/180);
        st = sin(deg*pi/180);
        new = [ (a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-ct)+x*ct+(-c*v+b*w-w*y+v*z)*st
        (b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-ct)+y*ct+(c*u-a*w+w*x-u*z)*st
        (c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-ct)+z*ct+(-b*u+a*v-v*x+u*y)*st ]';
        new = new/norm(new);
        %quiver3(a,b,c,new(1),new(2),new(3))

%         syms tt
        limit2 = 10;%solve(sqrt(sum((new*tt).^2))==3);

        for tt = 0:1:limit2 % inch along the line
            edgepoint =[a b c] + new*tt;
            x_new = edgepoint(1);
            y_new = edgepoint(2);
            z_new = edgepoint(3);

            i = round(x_new);%1 + (x_new - srow_x(4))/srow_x(3)
            j = round(y_new);%1 + (y_new - srow_y(4))/srow_y(1)
            k = round(z_new);%1 + (z_new - srow_z(4))/srow_z(2)

            % Search adjacent voxels - 9x3
            adjacent = [i-1     j-1     k
                        i-1     j       k
                        i-1     j+1     k
                        i       j-1     k
                        i       j       k
                        i       j+1     k
                        i+1     j-1     k
                        i+1     j       k
                        i+1     j+1     k];
            side_ones = [zeros(9,1) zeros(9,1) ones(9,1)];
            adjacent_plusk = adjacent + side_ones;% slice above
            adjacent_minusk = adjacent - side_ones;% slice below
            total_adjacent = [adjacent_plusk; adjacent; adjacent_minusk];% 27x3 matrix

            i_adj = total_adjacent(:,1);
            j_adj = total_adjacent(:,2);
            k_adj = total_adjacent(:,3);

            x_adj = i_adj;%(srow_x(1).*i_adj) .+ (srow_x(2).*j_adj) + (srow_x(3).*(k_adj-1)) .+ srow_x(4); %27x1
            y_adj = j_adj;%(srow_y(1).*(i_adj-1)) .+ (srow_y(2).*j_adj) .+ (srow_y(3).*k_adj) .+ srow_y(4);
            z_adj = k_adj;%(srow_z(1).*i_adj) .+ (srow_z(2).*(j_adj-1)) .+ (srow_z(3).*k_adj) .+ srow_z(4);

            x_diff = x_new - x_adj; % 27x1
            y_diff = y_new - y_adj;
            z_diff = z_new - z_adj;
            xyz_adj_dist = sqrt(x_diff.^2+y_diff.^2+z_diff.^2); %27x1
            pixrez=1.5;
            idx = find(xyz_adj_dist < pixrez);
            for rr = idx
                mask(i_adj(rr),j_adj(rr),k_adj(rr))=1; % for close voxels, mask = 1
            end
        end
    end

end


for pp = 72:88
   figure
   imagesc(squeeze(mask(:,pp,:)))
   axis off
   colormap gray
end
close all

% Save as NII.gz file
% "Close mask" using ITK afterwards









file = load_nii([script01_prefix '07_HistMatch/Brain_HistMatch_' T2FLAIR_SeriesInstanceUID '.nii.gz'])

dim = [file.hdr.dime.dim(4) file.hdr.dime.dim(3) file.hdr.dime.dim(2)]; % # pixels
origin_x = file.hdr.hist.qoffset_x % origin point - upper left
origin_y = file.hdr.hist.qoffset_y
origin_z = file.hdr.hist.qoffset_z
srow_x = file.hdr.hist.srow_x;
srow_y = file.hdr.hist.srow_y;
srow_z = file.hdr.hist.srow_z;



% 512x512x158 --> mm x mm x mm
i=1;
j=1;
k=1;
x = (srow_x(1)*i) + (srow_x(2)*j) + (srow_x(3)*(k-1)) + srow_x(4)
y = (srow_y(1)*(i-1)) + (srow_y(2)*j) + (srow_y(3)*k) + srow_y(4)
z = (srow_z(1)*i) + (srow_z(2)*(j-1)) + (srow_z(3)*k) + srow_z(4)
coords = [x y z]


% Get mm x mm x mm coordinates for center


% Get trajectory line based off 2 points
point1 = [-15.03 -53.67 58.24];
point2 = [-14.37 -52.98 46.27];
point3 = [-13.58 -52.37 40.3];

vector1 = point1-point2;
vector2 = point2-point3;
vector3 = point1-point3;

point2 = [1 2 3];
point1 = [0 0 0];

 
% Get VOI +/- 5 mm coords
vector_line = point2-point1;
syms t
sol = vpasolve(sqrt(sum(((vector_line*t).^2)))==5) % t | dist = 5
% dist = sqrt(sum(((vector*sol).^2))) ; % Check dist = 5

% End Points
xyz1 = point1 + vector_line*sol; % Get coords for which dist = 5 = h/2
xyz2 = point1 + vector_line*-sol; 

% Find circumference perpendicular plan
% Perp vector = [1 1 -(A+B)/C]
vector_line_unit = vector_line/norm(vector_line);
boundary_point = xyz1 + vector_perp1; % Rotate norm







