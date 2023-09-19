function out_pos = equiangular_coordinates(ele_min, ele_max, min_azi, max_azi, ele_res, azi_re)

% ele_min=-90; ele_max=90; min_azi=0; max_azi=360; ele_res=5; azi_re=5;
radius = 1; 
base_azi = min_azi:azi_re:max_azi-azi_re;
base_ele = sort([0:ele_res:ele_max, 0-ele_res:-ele_res:ele_min]);

ele = repmat(base_ele, [1, length(base_azi)]);
azi = sort(repmat(base_azi, [1, length(base_ele)]));
r = ones(length(ele(:)), 1) * radius;
raw_pos = [azi(:), ele(:), r];
raw_pos = unique(sortrows(raw_pos, 1), 'rows');

% ensure cartesian coordinates are unique
[x, y, z] = sph2cart(deg2rad(raw_pos(:, 1)), ...
                     deg2rad(raw_pos(:, 2)), ...
                     raw_pos(:, 3));
[~, idx] = unique([round(x, 2), round(y, 2), round(z, 2)], 'rows');
out_pos = raw_pos(idx, :);
out_pos = sortrows(out_pos,2);
end
% 

% scatter(o ut_pos(:, 1), out_pos(:, 2))