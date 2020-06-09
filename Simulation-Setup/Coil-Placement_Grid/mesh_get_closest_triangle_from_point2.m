function idx=mesh_get_closest_triangle_from_point2(m, point, varargin)
  % Gives the linear index to the closest triangle of the surface specified by triangle_regions to the point.
  % USAGE: idx=MESH_GET_CLOSEST_TRIANGLE_FROM_POINT(m, point[, triangle_regions=1])
  % Mirko Windhoff, 2010
  % Modified by Katie Mantell, 2020
  % $Id: mesh_get_closest_triangle_from_point.m 671 2011-05-06 16:22:49Z thielscher $
  
  if nargin>2, triangle_region=varargin{1}; else triangle_region=1005; end;
  m = mesh_extract_regions(m,'elemtype','tri','region_idx',triangle_region);
  center=(m.nodes(m.triangles(:,1),:) + ...
      m.nodes(m.triangles(:,2),:) + ...
      m.nodes(m.triangles(:,3),:))./3;
  [minval,idx]=min(sqrt(sum((center-repmat(point, size(center,1), 1)).^2,2)));
end
