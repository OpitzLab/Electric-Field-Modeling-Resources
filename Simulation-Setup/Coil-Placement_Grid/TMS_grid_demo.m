%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Demo for building a 3x3 TMS grid for SimNIBS simulations
% 
% Presented at the online workshop:
% Computational Modeling in Non-Invasive Brain Stimulation (NIBS) Workshop
% May 27-28, 2020
% 
% Katie Mantell, 2020
% Opitz Lab, University of Minnesota
% 
% for simnibs version 3.1.2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc
%% Data to enter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% ENTER/MODIFY DATA IN THIS SECTION %%%%%%%%%%%%%%%%%%%%

% these are inputs that I have found to be useful in setting up simulations
% for both hemispheres for multiple participants. you should modify these
% to work with your naming conventions and folder structures.
sim_directory = 'C:/Users/k/Desktop/Workshop_Demo/Simulations';         % parent folder for simulation folders
mesh_path = 'C:/Users/k/Desktop/Workshop_Demo/simnibs_examples/ernie';  % folder where head mesh model is
ID = 'ernie';                                                           % participant ID
date = '20200605';
hemisphere = 'left';

% head coordinates space: position of middle of grid (point in GM, skin, 
% etc, all ok). select using gmsh, SimNIBS GUI, EEG cap file, etc.
centerpoint_init = [-67.158218, -16.599871, 82.762589]; % ernie C3

% set up TMS grid
Nxpoints = 3;       % number of grid points in x direction
Nypoints = 3;       % number of grid points in y direction
spacing = 10;       % grid point spacing (in mm)
coil_theta = 45;	% rotation for coil at each position (in degrees)

% setting up direction vectors for different coordinate spaces.
% In most cases these should work. Most likely the only change you would
% want is the grid angled on the head, which is controlled by gridx
gridx = [1 0 0];    % grid coordinate space: choose direction for grid x axis ([1 0 0] for straight grid)
gridz = [0 0 1];    % grid coordinate space: choose direction for grid z axis (typically [0 0 1])
coilz = [0 0 1];    % coil coordinate space: choose coil z axis (typically [0 0 1])
headz = [0 0 1];    % head coordinate space: choose head z axis (will be [0 0 1] for meshes created with headreco)

% location of a template file. I set up a TMS simulation with the GUI and 
% ran it. I then kept the generated .mat file as a template for setting up
% future simulations
sim_template = 'C:/Users/k/Desktop/Workshop_Demo/simulationTemplate.mat'; 

% coil files come with SimNIBS. The .ccd files are nice for visualizing the
% coil in gmsh, however the .nii.gz files are faster for simulations
coil_file = 'C:/Users/k/AppData/Local/SimNIBS/miniconda3/envs/simnibs_env/Lib/site-packages/simnibs/ccd-files/Magstim_70mm_Fig8.ccd';

%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ENTRY COMPLETE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create paths and file name
% modify these as needed for your naming conventions and folder structures
sim_name = sprintf('%s_%s_%s',ID,date,hemisphere); % simulation folder name
sim_folder = sprintf('%s/%s',sim_directory,sim_name); % simulation folder path

% create new simulation folder if it does not exist already
% NOTE: you have to be able to run matlab as an administrator in Windows
% alternative is to comment out these 3 lines and create the folder by hand
if ~exist(sim_folder,'dir')
    mkdir(sim_folder)
end

%% create TMS grid
% load mesh
m = mesh_load_gmsh4(sprintf('%s/%s.msh',mesh_path,ID));

% region number for skin
skinRegion = 1005;

% get skin from mesh
skin = mesh_extract_regions(m,'elemtype','tri','region_idx',skinRegion);

% get skin surface data (center points and normal vectors)
% centers are used because there is a normal vector associated with that
% point, however mesh nodes can also be used.
centers = mesh_get_triangle_centers(skin);
normals = mesh_get_triangle_normals2(skin);

% get data for skin point closest to chosen center point
idx=mesh_get_closest_triangle_from_point2(skin, centerpoint_init, skinRegion);
centerpoint = centers(idx,:);
centernormal = normals(idx,:);

% rotate grid axis based on center point normal vector
v = cross(gridz,centernormal);
theta = acos(dot(gridz,centernormal)/(norm(gridz)*norm(centernormal)));
c = dot(gridz,centernormal)*cos(theta);
v_skew = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R = eye(3) + v_skew + v_skew^2*(1/(1+c));
xaxis_new = (R*gridx')';
xaxis_new = xaxis_new/normd(xaxis_new);     % normalize vector
yaxis_new = cross(centernormal,xaxis_new);
yaxis_new = yaxis_new/normd(yaxis_new);     % normalize vector


% get grid positions and normal vectors
% this calculates a Ny x Nx grid around a center point on the scalp
points = zeros(Nxpoints,Nypoints,3);    % initialize matrix for coordinates
zvector = zeros(Nxpoints,Nypoints,3);   % initialize matrix for coil zvectors

% loop through grid using spacing indexing
for i = -floor(Nxpoints/2):ceil(Nxpoints/2)-1
    for j = -floor(Nypoints/2):ceil(Nypoints/2)-1
        % get index of triangle closest to point estimated by spacing
        idx = mesh_get_closest_triangle_from_point2(skin,(centerpoint + spacing*i * xaxis_new + spacing*j * yaxis_new));
        % save coordinates and normal vector in matrix with matrix indexing
        points(i+floor(Nxpoints/2)+1,j+floor(Nypoints/2)+1,:) = centers(idx,:);
        zvector(i+floor(Nxpoints/2)+1,j+floor(Nypoints/2)+1,:) = normals(idx,:);
    end
end

% go from Ny x Nx x 3 to Number of positions x 3 ([x1,y1,z1;x2,y2,z2;...]) matrices
grid = reshape(points,Nxpoints * Nypoints,3);
zvector = reshape(zvector,Nxpoints * Nypoints,3);

%% get coil y vector

% calculating coil yvector
% 1st: project headz onto coil xy-plane (zvector is normal to this plane)
coily = repmat(headz,length(zvector),1) - (dot(repmat(headz,length(zvector),1),zvector,2)./dot(zvector,zvector,2)).*zvector;

yvector = zeros(size(coily));

% 2nd: rotate projected vector by theta
for i = 1:length(coily)
    ux = zvector(i,1);
    uy = zvector(i,2);
    uz = zvector(i,3);
    R = [cos(theta)+ux^2*(1-cos(theta)), ux*uy*(1-cos(theta))-uz*sin(theta), ux*uz*(1-cos(theta))+uy*sin(theta);...
        uy*ux*(1-cos(theta))+uz*sin(theta), cos(theta)+uy^2*(1-cos(theta)), uy*uz*(1-cos(theta))-ux*sin(theta);...
        uz*ux*(1-cos(theta))-uy*sin(theta), uz*uy*(1-cos(theta))+ux*sin(theta), cos(theta)+uz^2*(1-cos(theta))];
    yvector(i,:) = (R*coily(i,:)')';
end

%% normalize vectors

% normalize coil y vectors
yvector = bsxfun(@rdivide,yvector,normd(yvector));

% get normalized coil x vector from z and y vectors
xvector = cross(zvector, yvector);
xvector = bsxfun(@rdivide,xvector,normd(xvector));


%% update simulation mat file template

% set number of simulations
nTMS = length(yvector);

% load template simulation .mat file
t=load(sim_template);

% get template TMS poslist
TMS = t.poslist{1}.pos(1);

% set head mesh data
% see simnibs website for more info on simulation structures:
% https://simnibs.github.io/simnibs/build/html/documentation/sim_struct/sim_struct.html
t.eeg_cap = sprintf('%s/m2m_%s/eeg_positions/EEG10-10_UI_Jurak_2007.csv',mesh_path,ID); % location of EEG cap from headreco
t.fnamehead = sprintf('%s/%s.msh',mesh_path,ID);                                        % location of head mesh
t.pathfem = sim_folder;                                                                 % simulation folder path
t.fields = 'eEjJ';                                                                      % which fields should be saved
t.subpath = sprintf('%s/m2m_%s',mesh_path,ID);                                          % m2m folder path (this folder is generated in headreco)
t.fname_tensor = '';                                                                    % location of tensor file if one exists
t.open_in_gmsh = 0;                                                                     % do not automatically open gmsh at end of simulation
t.poslist{1}.fnamecoil = coil_file;                                                     % location of coil file

% set each position in poslist by looping through grid points
for i = 1:nTMS
    % add position i to poslist using TMS template
    t.poslist{1}.pos(i) = TMS;
    % write in coil vector and position (coil is moved out 4mm in z
    % direction to account for hair,plastic around coil, etc.)
    t.poslist{1}.pos(i).matsimnibs = cat(2,[xvector(i,:)'; 0],[yvector(i,:)'; 0],[-zvector(i,:)'; 0],[grid(i,:)' + (4 * zvector(i,:))'; 1]);
end

%% save simulation mat file and data for visualization

% the saved .mat file can be run in the command line:
% simnibs <simulation file>.mat
save(sprintf('%s/%s.mat',sim_folder,sim_name),'-struct','t');

% save the coil grid positions and coil vectors. Can be useful for later
% analyses
save(sprintf('%s/coil_data.mat',sim_folder),'grid','xvector','yvector','zvector');

% visualization in gmsh is a useful way to make sure grid and vectors look
% correct
write_coilorientation_in_mesh(zvector,grid,sprintf('%s/coil_grid_%s_z.geo',sim_folder,hemisphere));
write_coilorientation_in_mesh(xvector,grid,sprintf('%s/coil_grid_%s_x.geo',sim_folder,hemisphere));
write_coilorientation_in_mesh(yvector,grid,sprintf('%s/coil_grid_%s_y.geo',sim_folder,hemisphere));
write_pointvalues_in_geofile(grid,1:nTMS,sprintf('%s/coil_pos_nums_%s',sim_folder,hemisphere));