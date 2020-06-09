function write_coilorientation_in_mesh(vec, points, viewfilename)
  % writes the coil orientation vector(s) in the mesh.
  % USAGE:
  % WRITE_VECTOR_IN_MESH(vec, points, viewfilename)
  % where size(points)=[n 3]
  % Mirko Windhoff, 2009
  % $Id: points_write_V1_view.m 481 2010-09-15 15:50:21Z aopitz $
  tic
  fd=fopen(viewfilename,'wt');
  fwrite(fd, ['View "' viewfilename '" {']);
  for i=1:size(points,1)
    fwrite(fd, sprintf('VP (%f, %f, %f) {%f, %f, %f};\n', points(i,:), vec(i,:)));
  end;
  fwrite(fd, '};');
  fclose(fd);
  toc
end