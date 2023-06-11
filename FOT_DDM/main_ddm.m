%% Data generator for code validation wrt Monte Carlo
solver_options.itersolver = 'none'; 
% solver_options.restart = 30;
% solver_options.tol = 1e-12; 
% solver_options.maxit=50;


addpath('problemsetup\')
addpath('forwardsolver\')
addpath('smw\');



assemble_on_fly=1;
SettingsFile = 'test_el_3D_flo';
%%
disp('Running for basic SPN with source at 1mfp')
func = @FwdSolver;

run(SettingsFile) %Contains the Problem setting

NofSources = size(Source,2); % Here Source is a structure-type variable

% AOI is specified as angle wrt normal outside the medium. Henceforth AOI
% refers to the angle inside the medium and wrt to the outward normal at
% the corresponding interface.
for s=1:NofSources
    if Source(s).AOI ~=90
        Source(s).AOI =pi - asin(Source(s).nfibre*sin(Source(s).AOI*pi/180)/Background_property(1,end));
        Source(s).AOI=Source(s).AOI*180/pi; % Store angle in degrees
    end
    Source(s).mfreq = 2*pi*Source(s).mfreq*1e6; % Store modulation frequency
end

%%

%%%%% Commented on 26th June 2022 %%%%%%%%%%%%%%
if (Dimension==2)
    if isempty(shape)
        [mesh,Src,Det]=MeshGen2D_new(saveMesh,Domain,h,Source,Detector); % Generate the triangulation and meshing information
    else
        [mesh,Src,Det]=MeshGen2D_new(saveMesh,Domain,h,Source,Detector,shape,centre,r);
    end
    
    % This part works only for homogeneous meshes. The Assembly codes will have
    % to be correspondingly modified for inhomogeneous/ adaptive meshes
    if ~assemble_on_fly
        AssemblyFile=Assemble2D_new(mesh);
    else
        AssemblyFile = @Assemble_on_fly;
    end
else
    if isempty(shape)
        [mesh,Src,Det]=MeshGen3D_new(saveMesh,Domain,h,Source,Detector); % Generate the triangulation and meshing information
    else
        [mesh,Src,Det]=MeshGen3D_new(saveMesh,Domain,h,Source,Detector,shape,centre,r);
    end
    
    % This part works only for homogeneous meshes. The Assembly codes will have
    % to be correspondingly modified for inhomogeneous/ adaptive meshes
    if ~assemble_on_fly
        AssemblyFile=Assemble3D(mesh);
    else
        AssemblyFile = @Assemble_on_fly_3D;
    end
end
%%

%  mesh=MeshGen3D_new(Domain,h);
 
NofElem=size(mesh.tri,1); % # of elements
n0=external_ref_ind;

%mesh.opt is sub-structure that contains the optical properties of the
%medium
if simulation =='F'
    % For fluorescence simulations
    
    mesh.opt.muaxi=Background_property(1)*ones(1,1,NofElem);%intrinsic absorption coeffecient at excitation wavelength
    mesh.opt.musx=Background_property(2)*ones(1,1,NofElem); % scattering coefficient at excitation wavelength
    mesh.opt.muaimult=Background_property(3); % multiplication factor for intrinsic absorption coefficient at emission wavelength
    mesh.opt.musm=Background_property(4)*mesh.opt.musx; % scattering coefficient at emission wavelength
    mesh.opt.gx=Background_property(5); % anisotropy factor at excitation
    mesh.opt.gm=Background_property(6); %anisotropy factor at emission
    mesh.opt.refrind=Background_property(7); % refractive index of the medium
    mesh.opt.speed=speed_of_light/mesh.opt.refrind; % speed of light in the medium
    mesh.opt.refrind_out =n0; %external refractive index
    
    mesh.opt.gamma = gammaval;  % gamma parameter for modified HG
    mesh.opt.alpha = alphaval; % alpha parameter for TTHG
    mesh.opt.gbx = gb_x; % backscattering anisotropy fot TTHG
    mesh.opt.gbm = gb_m; % backscattering anisotropy fot TTHG
    %Data for fluorophore
    %     muaxf=zeros(1,1,NofElem);
    mesh.opt.muaxf=Fluorophore_property(1)*ones(1,1,NofElem); % fluorophore absorption coefficient at excitation wavelength
    % if there is a inclusion, assign muaxf values to it
    if ~isempty(shape)
        index_inhom2 = [0,cumsum(mesh.inhom.ind)];
        
        for ii=1:length(Inhom_property)
            %     mesh.opt.muaxf(1,1,tF(1+index_inhom2(ii)))= Inhom_property(1);
            mesh.opt.muaxf(1,1,mesh.inhom.elem(1+index_inhom2(ii):index_inhom2(ii+1)))=Inhom_property(ii);
        end
    end
    mesh.opt.muafmult=Fluorophore_property(2); %  multiplication factor for fl abs coeff at emission wavelength
    mesh.opt.Qf=Fluorophore_property(4)*ones(1,1,NofElem); % quantum efficiency of the fluorophore
    mesh.opt.tau=Fluorophore_property(5); % fluorescence lifetime in ns
    
    mesh.opt.bngx = mesh.opt.gx;
    mesh.opt.bngm = mesh.opt.gm;
    SimType=1; % Set to 1 for fluorescence simulation, else set to 0
    %     Plot_mu(mesh.nodes,mesh.tri,mesh.hx,mesh.opt.muaxf);
    
else
    %Data for Layer 1
    mesh.opt.muaxi=Background_property(1)*ones(1,1,NofElem); % absorption coefficient
    mesh.opt.musx=Background_property(2)*ones(1,1,NofElem); %scattering coefficient
    % if there is an inclusion, assign muaxi values to it
    if ~isempty(shape)
        index_inhom2 = [0,cumsum(mesh.inhom.ind)];
        for ii=1:length(Inhom_property)
            mesh.opt.muaxi(1,1,mesh.inhom.elem(1+index_inhom2(ii):index_inhom2(ii+1)))=Inhom_property(ii);
        end
    end  
    
    mesh.opt.gx=Background_property(3); %anisotropy factor
    mesh.opt.refrind=Background_property(4); % refractive index of the medium
    mesh.opt.speed=speed_of_light/mesh.opt.refrind; % speed of light in medium
    mesh.opt.refrind_out =n0; % speed of light outside
    
    mesh.opt.gamma = gammaval;  % gamma parameter for modified HG
    mesh.opt.alpha = alphaval; % alpha parameter for TTHG
    mesh.opt.gbx = gb_x; % backscattering anisotropy fot TTHG
  
    mesh.opt.bngx = mesh.opt.gx;
    SimType=0;% Set to 1 for fluorescence simulation, else set to 0
    %     Plot_mu(mesh.nodes,mesh.tri,mesh.hx,mesh.opt.muaxi);
end


% The angle of incidence is specified wrt the outward normal on a given
% interface. For the computations, we need to convert this to angle wrt
% global coordinate system that is wrt y/z axis. 

for s=1:NofSources
    if strcmp(Src(s).Type,'collimated')
        temp =180-Src(s).AOI;
        if Src(s).Loc(1,1)==Domain(1)
            Src(s).ang_in = 90+temp;
        elseif Src(s).Loc(1,1)==Domain(2)
            Src(s).ang_in = 270+temp;
        elseif Src(s).Loc(1,2)==Domain(3)
            Src(s).ang_in =temp; 
        elseif Src(s).Loc(1,2)==Domain(4)
            Src(s).ang_in=360-Src(s).AOI;
        elseif Src(s).Loc(3)==Domain(5)
                %         normal_in = '+z'
                Src(s).ang_in =temp; %180-Src(s).AOI ;
        elseif Src(s).Loc(3)==Domain(6)
                Src(s).ang_in=360-Src(s).AOI;
        end
            %mu_t_mean = mean(mesh.opt.muaxi + mesh.opt.muaxf + mesh.opt.musx*(1-mesh.opt.gx));
            mu_t_mean = mean(mesh.opt.muaxi + mesh.opt.musx*(1-mesh.opt.gx));
            mfp = 1/(mu_t_mean);
            Src(s).ID = 'I';
            %           Src(s) = get_one_mfp_source(Domain,mesh,Src(s),mfp);
            Src(s) = get_one_mfp_source(Dimension,mesh,Src(s),mfp);
    end
end

% [Jd,PhiX]=func(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,pfuncflag,evalFluence,probType);

%%
mesh2 = mesh;   

solver_options.itersolver = 'none'; 
solver_options.restart = 30;
solver_options.tol = 1e-12; 
solver_options.maxit=50;
solver_options.pfuncflag = pfuncflag; 
solver_options.evalFluence = evalFluence; 
solver_options.probType = probType;


mesh_spacing=mean([mesh.hx,mesh.hy,mesh.hz]);

% alpha=(1/(0.25*mesh_spacing));
% beta=1/(2*(mesh_spacing^2));

alpha=0.1307;
beta=1;

NoDomains=2;
partmat= ddm_interfacing(mesh,NoDomains);

[Jd_pert,PhiX_pert]=FwdSolver_DDM(Dimension,N,SimType,mesh2,partmat,Src,Det,alpha,beta,AssemblyFile,solver_options);

