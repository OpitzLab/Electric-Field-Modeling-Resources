function [triangle_normals] = mesh_get_triangle_normals2(m)
  % Calcuates the normals for each triangles. If the orientation of the
  % triangle edges is ccw to a viewer, then the normal points to the viewer.
  % USAGE: triangle_normals=MESH_GET_TRIANGLE_NORMALS(m)
  % Mirko Windhoff, 2009
  % Modified by Ivan Alekseichuk, 2018
  
  a                = m.nodes(m.triangles(:,1),:) - m.nodes(m.triangles(:,2),:);
  b                = m.nodes(m.triangles(:,1),:) - m.nodes(m.triangles(:,3),:);
  triangle_normals = cross(a,b,2);
  triangle_normals = bsxfun(@rdivide, triangle_normals, normd(triangle_normals));
  
end
