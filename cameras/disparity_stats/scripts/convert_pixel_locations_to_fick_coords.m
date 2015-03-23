function [] = convert_pixel_locations_to_fick_coords(sensor_width,sensor_center,focal_length)

%elevation is major circles, neg = down, pos = up
elevation = -((1:sensor_width)-sensor_center).*(2*atand(.5/focal_length));

%azimuth in minor circles, neg = left, pos = right
azimuth = -((1:sensor_width)-sensor_center).*(2*atand(.5/focal_length));
keyboard

