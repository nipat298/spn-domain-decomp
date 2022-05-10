% Setting file 
%Problem Setting to be used for the simulation
speed_of_light = 3e10; % Sets scale of prob. presently cm/s , set appropriately for mm or other scalings
N=3; % Order of approximation
simulation ='E'; % 'F' for fluorescence or 'E' for elastic scattering
probType ='fwd'; % fwd- forward solution, 'recon' - store data for inverse algo
evalFluence =1; % Set to 1 to save Fluence values as well, else set to 0
saveMesh = 1; % Set to 1 to save the mesh generated else set to 0.
pfuncflag = 0; % type of phase function. set to 0 - Henyey Greenstein, 
%                                        1 - Modified Henyey Greenstein, 
%                                        2- two term Henyey Greenstein Phase function
%% Problem Geometry
Dimension = 2; % Choose '2' for 2D and '3' for 3D
h=[0.025,0.025]; %Mesh spacing for structured mesh. Recommended to use same mesh size in each dimension, 
                % however you can specify different mesh size if needed.
                % [hx,hy] or [hx,hy,hz]
Domain=[-1,1,-1,1]; %,0,0.1]; %[xmin,xmax,ymin,ymax] or [xmin,xmax,ymin,ymax,zmin,zmax] respectively
%Specifications of the inhomogeneities
shape = 'circle';
centre=[0.5,0];
   r=[0.25];

%% Medium properties
%Properties at the excitation wavelength
%Use property value 'NAN' if not applicable
external_ref_ind=[1]; %refractive index outside. NO mismatch. Probe inserted into the medium
%external_ref_ind at [top,right,bottom,left] interface; If only one val is
%satisfied, same val is assumed on all interfaces.
Background_property=[0.031,54.75,0.8,1.37];%[mua_i, mus,muai_mult, mus_mult,g_x,g_m,refr_ind]; % One row for each layer
Inhom_property=[0.3]; % Specify one row for each inhomogeneity
gammaval = 1.8; % gamma parameter for modified henyey greenstein phase function
alphaval =0.9; gb_x = 0.3; gb_m=0.3 ; % parameters for TTHG
%% Source properties
s=1;
Source(s).Type='point'; % Type of source point or collimated
Source(s).mfreq = 100; % Modulation frequency
Source(s).AOI = 90; % angle of incidence on the boundary wrt normal. Set to 90 for Isotropic source
Source(s).ID='S'; % Choose 'S' for surface, 'I' for internal source
Source(s).Strength=1; % Source strength in mW
Source(s).nfibre =1; % refractive index of fibre core
Source(s).NA =0.22; % numerical aperture of the fibre
Source(s).Loc=[0,1]; % source location
Source(s).Radius = 0.000001; % set beam radius. For pencil beam set radius less than min(hx,hy)

s=2;
Source(s).Type='point'; % Type of source point or collimated
Source(s).mfreq = 100; % Modulation frequency
Source(s).AOI = 90; % angle of incidence on the boundary wrt normal. Set to 90 for Isotropic source
Source(s).ID='S'; % Choose 'S' for surface, 'I' for internal source
Source(s).Strength=1; % Source strength in mW
Source(s).nfibre =1; % refractive index of fibre core
Source(s).NA =0.22; % numerical aperture of the fibre
Source(s).Loc=[1,0]; % source location
Source(s).Radius = 0.000001; % set beam radius. For pencil beam set radius less than min(hx,hy)

s=3;
Source(s).Type='point'; % Type of source point or collimated
Source(s).mfreq = 100; % Modulation frequency
Source(s).AOI = 90; % angle of incidence on the boundary wrt normal. Set to 90 for Isotropic source
Source(s).ID='S'; % Choose 'S' for surface, 'I' for internal source
Source(s).Strength=1; % Source strength in mW
Source(s).nfibre =1; % refractive index of fibre core
Source(s).NA =0.22; % numerical aperture of the fibre
Source(s).Loc=[0,-1]; % source location
Source(s).Radius = 0.000001; % set beam radius. For pencil beam set radius less than min(hx,hy)

s=4;
Source(s).Type='point'; % Type of source point or collimated
Source(s).mfreq = 100; % Modulation frequency
Source(s).AOI = 90; % angle of incidence on the boundary wrt normal. Set to 90 for Isotropic source
Source(s).ID='S'; % Choose 'S' for surface, 'I' for internal source
Source(s).Strength=1; % Source strength in mW
Source(s).nfibre =1; % refractive index of fibre core
Source(s).NA =0.22; % numerical aperture of the fibre
Source(s).Loc=[-1,0]; % source location
Source(s).Radius = 0.000001; % set beam radius. For pencil beam set radius less than min(hx,hy)


%% Detector properties
Detector.AOD = 90;% angle of detection. set to 90 for isotropic detector
Detector.NA = 'NAN'; % numerical aperture of the detector fibre
Detector.nfibre =1; % refractive index of the detector fiber core
Detector.Type= 'point'; % TYpe of detector
Detector.ID='S'; % Choose 'S' for surface detector, 'I' for internal detectors

%% Probe Geometry
% ax=[Domain(1)+0.05:0.05:-0.05,0.05:0.05:Domain(2)-0.05]';
% ay=[Domain(3)+0.05:0.05:0.45,0.55:0.05:Domain(4)-0.05]';
% ay = ax; 
% ax=[Domain(1)+0.005:0.01:-0.03,0.03:0.01:Domain(2)-0.005]';
% Detector.Loc=[ax,Domain(3)*ones(size(ax));
%               Domain(1)*ones(size(ax)),ay;
%               ax,Domain(4)*ones(size(ax));
%               Domain(2)*ones(size(ax)),ay];
%               
ax=[-0.9:0.15:0.95]';
Detector.Loc=[ax,-1*ones(size(ax))
        ones(size(ax)),ax
        -ax,ones(size(ax))
        -1*ones(size(ax)),-ax];