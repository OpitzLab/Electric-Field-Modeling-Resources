function write_pointvalues_in_geofile(points,vals,fname)
% writes a Nx3 point list with Nx1 values in an .geo file which is named fname.geo
% USAGE: write_points_in_geofile(points,vals,fname)
% Modified by Katie Mantell, 2018

  
geofname = [fname '.geo'];
fd = fopen(geofname, 'w');
fwrite(fd, ['View "' geofname '" {']);
for j=1:length(vals)
    fwrite(fd, sprintf('SP (%f, %f, %f) {%f};\n', points(j,:), vals(j)));
end
fwrite(fd, '};');
fclose(fd);

end