% Setting file 
%Problem Setting to be used for the simulation
speed_of_light = 3e10; % Sets scale of prob. presently cm/s , set appropriately for mm or other scalings
N=3; % Order of approximation
simulation ='F'; % 'F' for fluorescence or 'E' for elastic scattering
probType ='recon'; % fwd- forward solution, 'recon' - store data for inverse algo
evalFluence =0; % Set to 1 to save Fluence values as well, else set to 0
saveMesh = 1; % Set to 1 to save the mesh generated else set to 0.
pfuncflag = 0; % type of phase function. set to 0 - Henyey Greenstein, 
%                                        1 - Modified Henyey Greenstein, 
%                                        2- two term Henyey Greenstein Phase function
%% Problem Geometry
Dimension = 3; % Choose '2' for 2D and '3' for 3D

h=[0.05,0.05,0.05]; %Mesh spacing for structured mesh. Recommended to use same mesh size in each dimension, 
                % however you can specify different mesh size if needed.
                % [hx,hy] or [hx,hy,hz]
Domain=[-1,1,-1,1,0,2];  %[xmin,xmax,ymin,ymax] or [xmin,xmax,ymin,ymax,zmin,zmax] respectively
%Specifications of the inhomogeneities
shape = [];

%% Medium properties
%Properties at the excitation wavelength
%Use property value 'NAN' if not applicable
external_ref_ind=[1,1.37,1.37,1.37,1.37,1.37]; %refractive index outside. NO mismatch. Probe inserted into the medium
%external_ref_ind at [top,right,bottom,left] interface; If only one val is
%satisfied, same val is assumed on all interfaces.
Background_property=[0.1,100,0.9,1.37];%[mua_i, mus,muai_mult, mus_mult,g_x,g_m,refr_ind]; % One row for each layer
Inhom_property=[0;0]; % Specify one row for each inhomogeneity
gammaval = 1.8; % gamma parameter for modified henyey greenstein phase function
alphaval =0.9; gb_x = 0.3; gb_m=0.3 ; % parameters for TTHG
%% Source properties
Source(1).Type='point'; % Type of source point or collimated
Source(1).mfreq = 100; % Modulation frequency
Source(1).AOI = 90; % angle of incidence on the boundary wrt normal. Set to 90 for Isotropic source
Source(1).ID='S'; % Choose 'S' for surface, 'I' for internal source
Source(1).Strength=complex(cos(0.15),sin(0.15)); % Source strength in mW
Source(1).nfibre =1; % refractive index of fibre core
Source(1).NA =0.22; % numerical aperture of the fibre
Source(1).Loc=[-0.3,0,0]; % source location
%Source(1).Loc=[-1,-1,0];
Source(1).Radius = 0.000001; % set beam radius. For pencil beam set radius less than min(hx,hy)
% 
% Source(2).Type='point'; % Type of source point or collimated
% Source(2).mfreq = 0; % Modulation frequency
% Source(2).AOI = 90; % angle of incidence on the boundary wrt normal. Set to 90 for Isotropic source
% Source(2).ID='S'; % Choose 'S' for surface, 'I' for internal source
% Source(2).Strength=complex(cos(0.15),sin(0.15)); % Source strength in mW
% Source(2).nfibre =1; % refractive index of fibre core
% Source(2).NA =0.22; % numerical aperture of the fibre
% Source(2).Loc=[0,0.25,0]; % source location
% Source(2).Radius = 0.000001; % set beam radius. For pencil beam set radius less than min(hx,hy)

%% Detector properties
Detector.AOD = 90;% angle of detection. set to 90 for isotropic detector
Detector.NA = 'NAN'; % numerical aperture of the detector fibre
Detector.nfibre =1; % refractive index of the detector fiber core
Detector.Type= 'point'; % TYpe of detector
Detector.ID='S'; % Choose 'S' for surface detector, 'I' for internal detectors

%% Probe Geometry
% ax=[Domain(1):h(1):-0.01,0.01:h(1):Domain(2)]';
ax=[Domain(1)+0.005:0.01:-0.01,0.01:0.01:Domain(2)-0.005]';

% ax=[Domain(1)+0.005:0.01:-0.03,0.03:0.01:Domain(2)-0.005]';
Detector.Loc=[ax,0*ones(size(ax)),0*ones(size(ax))];
              
