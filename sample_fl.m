% Setting file 
%Problem Setting to be used for the simulation
speed_of_light = 3e10; % Sets scale of prob. presently cm/s , set appropriately for mm or other scalings
N=3; % Order of approximation
simulation ='F'; % 'F' for fluorescence or 'E' for elastic scattering
probType ='fwd'; % fwd- forward solution, 'recon' - store data for inverse algo
evalFluence =1; % Set to 1 to save Fluence values as well, else set to 0
saveMesh = 1; % Set to 1 to save the mesh generated else set to 0.
pfuncflag = 0; % type of phase function. set to 0 - Henyey Greenstein, 
%                                        1 - Modified Henyey Greenstein, 
%                                        2- two term Henyey Greenstein Phase function
%% Problem Geometry
Dimension = 2; % Choose '2' for 2D and '3' for 3D
h=[0.01,0.01]; %Mesh spacing for structured mesh. Recommended to use same mesh size in each dimension, 
                % however you can specify different mesh size if needed.
                % [hx,hy] or [hx,hy,hz]
Domain=[-0.5,0.5,0,1]; %,0,0.1]; %[xmin,xmax,ymin,ymax] or [xmin,xmax,ymin,ymax,zmin,zmax] respectively
%Specifications of the inhomogeneities
% shape = [];
shape='circle'; % Type of inhomogeneity - 'circle', 'rect', 'ellipse', or 'bean' in 2D
                % 'sphere','cuboid' in 3D
% Specify centre of inhomogeneity. use one row for each inhomogeneity [xc,yc]
centre=[0,0.25
        0,0.75];
 % use this to specify 'r' = radius for circle
 % r=[l,b] for rect
 % r=[major axis, minor axis] for ellipse
 % r = scale for bean
 % Use one row for each inhomogeneity
r=[0.1
   0.1];

%% Medium properties
%Properties at the excitation wavelength
%Use property value 'NAN' if not applicable
external_ref_ind=[1,1,1.5,1.5]; %refractive index outside. NO mismatch. Probe inserted into the medium
%external_ref_ind at [top,right,bottom,left] interface; If only one val is
%satisfied, same val is assumed on all interfaces.
Background_property=[0.031,54.75,0.7987,0.732,0.8,0.8,1.37];%[mua_i, mus,muai_mult, mus_mult,g_x,g_m,refr_ind]; % One row for each layer
Fluorophore_property=[0.006,0.0846,1,0.016,0.56*1e-9]; %[mua_f, muaf_mult, ContrastRatio, Quant_eff, Fluor_lifetime]; %One row for each fluorophore 
Inhom_property=[0.3;0.3]; % Specify one row for each inhomogeneity
gammaval = 2; % gamma parameter for modified henyey greenstein phase function
alphaval =1; gb_x = 0; gb_m=0; % parameters for TTHG
%% Source properties
Source(1).Type='collimated'; % Type of source point or collimated
Source(1).mfreq = 100; % Modulation frequency
Source(1).AOI = 90; % angle of incidence on the boundary wrt normal. Set to 90 for Isotropic source
Source(1).ID='S'; % Choose 'S' for surface, 'I' for internal source
Source(1).Strength=1; % Source strength in mW
Source(1).nfibre =1; % refractive index of fibre core
Source(1).NA =0.22; % numerical aperture of the fibre
Source(1).Loc=[-0.25,0]; % source location

Source(2).Type='collimated';
Source(2).mfreq = 100;
Source(2).AOI = 0; % angle of incidence on the boundary
Source(2).ID='S'; % Choose 'S' for surface, 'I' for internal source
Source(2).Strength=1;
Source(2).nfibre =1;
Source(2).NA =0.22;
Source(2).Loc=[0.25,0];

%% Detector properties
Detector.AOD = 90;% angle of detection. set to 90 for isotropic detector
Detector.NA = 'NAN'; % numerical aperture of the detector fibre
Detector.nfibre =1; % refractive index of the detector fiber core
Detector.Type= 'point'; % TYpe of detector
Detector.ID='S'; % Choose 'S' for surface detector, 'I' for internal detectors

%% Probe Geometry
ax=[-0.45:0.2:0.45]';
ay=[0.05:0.2:0.95]';

Detector.Loc=[ax,0*ones(size(ax));
              0.5*ones(size(ay)),ay;
              ax,ones(size(ax));
              -0.5*ones(size(ay)),ay;];
       
