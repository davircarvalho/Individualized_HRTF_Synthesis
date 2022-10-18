function out_pos = equiangular_coordinates(ele_min, ele_max, min_azi, max_azi, ele_res, azi_re)
radius = 1; 
base_azi = min_azi:azi_re:max_azi-azi_re;
base_ele = sort([0:ele_res:ele_max , 0-ele_res:-ele_res:ele_min]);

ele = repmat(base_ele, [length(base_azi), 1]);
azi = repmat(base_azi, [1,length(base_ele)]);
r = ones(length(ele(:)), 1) * radius;
out_pos = [azi(:), ele(:), r];
end
