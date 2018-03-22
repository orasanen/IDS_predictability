function y = resample_by_interpolation(x,target_length)   

len = length(x);

target_points = linspace(1,len,target_length);
original_points = (1:len)';

y = interp1(original_points,x,target_points);
