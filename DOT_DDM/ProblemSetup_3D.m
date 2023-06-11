% DEFINE CONSTANTS
speed_of_light = 3*1e10; % speed of light in cm/s since all our units are in cms
mua = 0.031;
mus = 54.75; % Scattering coefficient
g = 0.8; % ANisotropy factor
n_in = 1.37; % refractive index of the medium 
n_out = 1; % refractive index of the outside medium
speed_in_med=speed_of_light/n_in; 
Ref_coeff_top = 0.431; % Reflection Coefficient at the top of the domain
A=(1+Ref_coeff_top)/(1-Ref_coeff_top);

mod_freq = 2*pi*100*1e6; % modulation frequency of the source in MHz. 

mfp = mean(1/(mua+(1-g)*mus)); % This is the transport mean free path

%% Source configuration

Source_Strength = 1;
% consider 2 sources on the boundary one at (0,-1) and other at (-1,0)
% Source_loc = [0,-1 + mfp;
%               -1+mfp,0
%               0, 1-mfp
%               1-mfp,0]; % (x,y) coordinates of the Source location. One row for each source 
          
 

% Source_loc = [0,-1 + mfp,0];

%Source_loc = [-1,1 ,0+mfp];
Source_loc = [-0.15,0,0];


          
% Detector configuration
% det_line(:,1) = -0.95:0.1:0.95; % One each edge the detectors are located at these points 
% % In most cases we do not want to place a detector at the source position.
% % The following matrix contains (x,y) coordinates of all the detector
% % positions
% NofDet_line = length(det_line);
% Det_loc =  [det_line,-1*ones(NofDet_line,1);
%             ones(NofDet_line,1),det_line;
%             -det_line,ones(NofDet_line,1);
%             -1*ones(NofDet_line,1),det_line];
