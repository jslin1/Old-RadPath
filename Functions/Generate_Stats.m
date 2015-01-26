function [o_mean o_std o_min o_max o_range o_numel] = Generate_Stats(i_values)
i_values = double(i_values);
o_mean = mean(i_values(:),1);
o_std = std(i_values(:),0,1);
o_min = min(i_values(:),[],1);
o_max = max(i_values(:),[],1);
o_range = o_max - o_min;
o_numel = size(i_values(:),1);
end

