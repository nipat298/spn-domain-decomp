Domain=[-1,1,-1,1,0,2]; % [xmin,xmax,ymin,ymax,zmin,zmax]
h=[0.05,0.05,0.05];
N=3;

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
