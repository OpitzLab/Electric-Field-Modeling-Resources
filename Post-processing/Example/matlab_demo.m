%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Demo for post-processing SimNIBS results in Matlab
% 
% Presented at the online workshop:
% Computational Modeling in Non-Invasive Brain Stimulation (NIBS) Workshop
% May 27-28, 2020
% 
% Sina Shirinpour, 2020
% Opitz Lab, University of Minnesota
% 
% for simnibs version 3.1.2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load mesh
% Load the mesh file by location
mesh = mesh_load_gmsh4('ernie_TMS_1-0001_Magstim_70mm_Fig8_scalar.msh');
%% Tissues indices
% 1: White Matter (WM) volume, 1001: White Matter surface
% 2: Grey Matter (GM) volume, 1002: Grey Matter surface
% 3: CSF volume, 1003: CSF surface
% 4: Bone volume, 1004: Bone surface
% 5: Skin volume, 1005: Skin surface
% 6: Eye balls, volume 1006: Eye balls surface
%% Grey matter extraction
% Extact all the GM tetrahedra
gm = mesh_extract_regions(mesh,'elemtype','tet','region_idx',2);
% Save extracted mesh for visualization
mesh_save_gmsh4(gm,'gm.msh');
%% Extract normE
% Find which structure in the element data corresponds normE data
field_idx = get_field_idx(gm,'normE','element');
% If there is no structure for what you want e.g. J in here, you get an error
% field_idx = get_field_idx(gm,'J','element');
% Extract all the normE data in the grey matter volume
gm_normE = gm.element_data{field_idx}.tetdata;
%% stats
% Calculate some statistics from the normE in grey matter
hist(gm_normE);
hist(gm_normE,50);
max_E = max(gm_normE) % @ dI/dt = 1 A/us
max_E_99 = prctile(gm_normE,99) % 99 percentile
mean_E = mean(gm_normE)
median_E = median(gm_normE)
%% Stats affected volume
% Take the values for the affected area (areas with at least half of the maximum E)
gm_normE_AV = gm_normE(gm_normE > max_E_99/2);
mean_E_AV = mean(gm_normE_AV)
median_E_AV = median(gm_normE_AV)
%% Basic automatic results
summary = mesh_get_fieldpeaks_and_focality(mesh,'field_idx','normE');
mesh_get_simulation_result
%% Interpolation
% Interpolate the results at any random location
% [X Y Z] coordinates of each point in each row
coords = [-34 -20 90; 34 -20 90];
% Interpolate the values at the desired coordinates in the mesh
iternp_values = get_fields_at_coordinates(mesh, coords);
%% Field map difference (question from previous day)
% Difference in electric field maps within a subject for two different montages
mesh1 = mesh_load_gmsh4('montage1.msh');
mesh2 = mesh_load_gmsh4('montage2.msh');
% Make sure meshes are the same size. If not, extract the GM first.
% If you want to compare multiple subjects, convert to a common space
% (e.g. MNI) first.

% Copy the mesh
mesh_diff = mesh1;
% Set the normE values of the new mesh as the difference of the normE in
% mesh1 and mesh2
mesh_diff.element_data{2}.tridata = mesh1.element_data{2}.tridata - ...
                                        mesh2.element_data{2}.tridata;
mesh_diff.element_data{2}.tetdata = mesh1.element_data{2}.tetdata - ...
                                        mesh2.element_data{2}.tetdata;
% Remove unnecessary fields (v and E here)
mesh_diff.element_data(1) = [];
mesh_diff.node_data = [];
% Save mesh
mesh_save_gmsh4(mesh_diff,'montage_diff.msh');