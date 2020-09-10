% Setting file 
%Problem Setting to be used for the simulation
speed_of_light = 3e10; % Sets scale of prob. presently cm/s , set appropriately for mm or other scalings
N=3; % Order of approximation
simulation ='F'; % 'F' for fluorescence
probType ='recon'; % fwd- forward solution, 'recon' - store data for inverse algo
evalFluence =1; % Set to 1 to save Fluence values as well, else set to 0
saveMesh = 1; % Set to 1 to save the mesh generated
pfuncflag = 0;
%% Problem Geometry
Dimension = 2; % Choose '2' for 2D and '3' for 3D
h=[0.05,0.05]; %Mesh spacing for structured mesh
Domain=[-0.5,0.5,0,1]; %,0,0.1]; %[xmin,xmax,ymin,ymax] or [xmin,xmax,ymin,ymax,zmin,zmax] respectively
    
%Specifications of the inhomogeneity
shape = [];
% shape='circle';
% centre=[0.05,0
%         -0.05,0];
% r=[0.02
%    0.02];

%% Medium properties
%Properties at the excitation wavelength
%Use property value 'NAN' if not applicable
external_ref_ind=[1,1.37,1.37,1.37]; %refractive index outside. NO mismatch. Probe inserted into the medium
Background_property=[0.031,54.75,0.7987,0.732,0.8,0.8,1.37];%[mua_i, mus,muai_mult, mus_mult,g_x,g_m,refr_ind]; % One row for each layer
Fluorophore_property=[0.006,0.0846,1,0.016,0.56*1e-9]; %[mua_f, muaf_mult, ContrastRatio, Quant_eff, Fluor_lifetime]; %One row for each fluorophore 
Inhom_property=[0.3;0.3]; % Specify one row for each inhomogeneity
gammaval = 1.8; % gamma parameter for modified henyey greenstein phase function
alphaval =1; gb_x = 0; gb_m=0; % parameters for TTHG
%% Source properties
% Source data structure. Indexed with Source(i). for individual entries
Source(1).Type='collimated';
Source(1).mfreq = 100;
Source(1).AOI = 0; % angle of incidence on the boundary
Source(1).ID='S'; % Choose 'S' for surface, 'I' for internal source
Source(1).Strength=1; 
Source(1).nfibre =1;
Source(1).NA =0.22;
Source(1).Loc=[0,0];

%% Detector properties
Detector.AOD = 90;
Detector.NA = 'NAN';
Detector.nfibre =1;
Detector.Type= 'point';
Detector.ID='S';

%% Probe Geometry
ax=[-0.5:0.05:0.5]';
Detector.Loc=[ax,0*ones(size(ax))];
       
