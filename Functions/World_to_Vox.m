function [i j k] = World_to_Vox(x, y, z, srow_x, srow_y, srow_z)

% World 1 (x), Vox 1/2/3
if srow_x(1)~=0  
    i = double(round(1 + (-x - srow_x(4))/srow_x(1))); % Vox 1 = i
elseif srow_x(2)~=0 
    j = double(round(1 + (-x - srow_x(4))/srow_x(2))); % Vox 2 = j
elseif srow_x(3)~=0 
    k = double(round(1 + (-x - srow_x(4))/srow_x(3))); % Vox 3 = k
end

% World 2 (y), Vox 1/2/3
if srow_y(1)~=0 
    i = double(round(1 + (-y - srow_y(4))/srow_y(1))); % Vox 1 = i
elseif srow_y(2)~=0 
    j = double(round(1 + (-y - srow_y(4))/srow_y(2))); % Vox 2 = j
elseif srow_y(3)~=0 
    k = double(round(1 + (-y - srow_y(4))/srow_y(3))); % Vox 3 = k
end

% World 3(z), Vox 1/2/3
if srow_z(1)~=0 
    i = double(round(1 + (z - srow_z(4))/srow_z(1))); % Vox 1 = i
elseif srow_z(2)~=0 
    j = double(round(1 + (z - srow_z(4))/srow_z(2))); % Vox 2 = j
elseif srow_z(3)~=0 
    k = double(round(1 + (z - srow_z(4))/srow_z(3))); % Vox 3 = k
end 


end