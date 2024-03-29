function [mesh,Md,varargout]=Datagen(simulation,SettingsFile,func,snr,fcs_flag,norm_flag,solver_options)
%% Data generator
% Run forward solver and obtain fluence and or partial current measurements
%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%
% simulation - type of simulation - 'E' - elastic, 'F' - fluorescence
% SettingsFile - script containing problem setting see sample_el.m (elastic) or
% sample_fl.m (fluorescence)
% func - function handle to forward solver
% snr (dB)- desired signal to noise ratio of measurements
% fcs_flag - first collision source flag '0' - basic spn, '1'- first
% collision source is used (fcs/delta-eddington variants)

%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%
% mesh - mesh structure as well as  optical properties
% Md - exiting partial current measurements with snr as specified


%%
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

if Dimension == 2 % Choose '2' for 2D and '3' for 3D
    if isempty(shape)
        [mesh,Src,Det]=MeshGen2D_new(saveMesh,Domain,h,Source,Detector); % Generate the triangulation and meshing information
    else
        [mesh,Src,Det]=MeshGen2D_new(saveMesh,Domain,h,Source,Detector,shape,centre,r);
    end
    
else
    % Not yet set up for 3D
    if isempty(shape)
        [mesh,Src,Det]=MeshGen3D(Domain,h,Src,Det);
    else
        [mesh,Src,Det]=MeshGen3D(Domain,h,Src,Det,shape,centre,r);
    end
end

% This part works only for homogeneous meshes. The Assembly codes will have
% to be correspondingly modified for inhomogeneous/ adaptive meshes
if Dimension==3
    AssemblyFile=Assemble3D(mesh.nodes,mesh.tri);
else if Dimension ==2
        AssemblyFile=Assemble2D_new(mesh);
    end
end
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
    
    SimType=1; % Set to 1 for fluorescence simulation, else set to 0
    %     Plot_mu(mesh.nodes,mesh.tri,mesh.hx,mesh.opt.muaxf);
    
else
    %Data for Layer 1
    mesh.opt.muaxi=Background_property(1)*ones(1,1,NofElem); % absorption coefficient
    mesh.opt.musx=Background_property(2)*ones(1,1,NofElem); %scattering coefficient
    % if there is an inclusion, assign muaxi values to it
    index_inhom2 = [0,cumsum(mesh.inhom.ind)];
    for ii=1:length(Inhom_property)
        mesh.opt.muaxi(1,1,mesh.inhom.elem(1+index_inhom2(ii):index_inhom2(ii+1)))=Inhom_property(ii);
    end
    mesh.opt.gx=Background_property(3); %anisotropy factor
    mesh.opt.refrind=Background_property(4); % refractive index of the medium
    mesh.opt.speed=speed_of_light/mesh.opt.refrind; % speed of light in medium
    mesh.opt.refrind_out =n0; % speed of light outside
    
     mesh.opt.gamma = gammaval;  % gamma parameter for modified HG
    mesh.opt.alpha = alphaval; % alpha parameter for TTHG
    mesh.opt.gbx = gb_x; % backscattering anisotropy fot TTHG
    
    SimType=0;% Set to 1 for fluorescence simulation, else set to 0
    %     Plot_mu(mesh.nodes,mesh.tri,mesh.hx,mesh.opt.muaxi);
end
% 
% if fcs_flag==1
%     % The angle of incidence is specified wrt the outward normal on a given
% % interface. For the computations, we need to convert this to angle wrt
% % global coordinate system that is wrt y/z axis. For 3D azimutthal symmetry
% % is implied in the SPn approximation. Once we have this, we can trace the
% % path of the ray from incidence till it exits the medium and obtain
% % corresponding element and node numbers as well as strength of the first
% % collimated source.
% 
%     for s=1:NofSources
%         temp =180-Src(s).AOI;
%         if Src(s).Loc(1,1)==Domain(1)
%             %         normal_in = '+x';
%             %         Src(s).ang_in = 270-Src(s).AOI;
%             Src(s).ang_in = 90+temp;
%         elseif Src(s).Loc(1,1)==Domain(2)
%             %         normal_in = '-x';
%             Src(s).ang_in = 270+temp;%90-Src(s).AOI;
%         elseif Src(s).Loc(1,2)==Domain(3)
%             %         normal_in = '+z'
%             Src(s).ang_in =temp; %180-Src(s).AOI ;
%         elseif Src(s).Loc(1,2)==Domain(4)
%             Src(s).ang_in=360-Src(s).AOI;
%         end
%         
%         if Src(s).Radius<min(mesh.hx,mesh.hy)
%         [Src(s).ray_elem,Src(s).seg_elem,Src(s).lseg_elem] = get_raytrace_2D(Domain,mesh.tri,mesh.nodes,Src(s).Loc(1,:),Src(s).elem,Src(s).ang_in*pi/180);
%         else
%             Src(s).beam = get_trace_2D_beam(Domain,mesh.tri,mesh.nodes,Src(s).Loc,Src(s).elem,Src(s).ang_in*pi/180,Src(s).Radius);
%         end
%     end
% else
%     mu_t_mean = mean(mesh.opt.muaxi + mesh.opt.muaxf + mesh.opt.musx*(1-mesh.opt.gx));
%     mfp = 1/(mu_t_mean);
%     % Adjust inocming angle and move sources 1mfp inside.
%     for s=1:NofSources
%         temp =180-Src(s).AOI;
%         
%         if Src(s).Loc(1,1)==Domain(1)
%             %         normal_in = '+x';
%             %         Src(s).ang_in = 270-Src(s).AOI;
%             Src(s).ang_in = 90+temp;
%             Src(s).Loc = Src(s).Loc + [mfp*sind(Src(s).ang_in),mfp*cosd(Src(s).ang_in)];
%             Src(s).ID = 'I';
%             Src(s).AOI = 90;
%             Src(s).elem = pointLocation(mesh.tri,Src(s).Loc(1),Src(s).Loc(2));
%              Src(s).node = nearestNeighbor(mesh.tri,Src(s).Loc(:,1),Src(s).Loc(:,2));
%         elseif Src(s).Loc(1,1)==Domain(2)
%             %         normal_in = '-x';
%             Src(s).ang_in = 270+temp;%90-Src(s).AOI;
%             Src(s).Loc = Src(s).Loc + [mfp*sind(Src(s).ang_in),mfp*cosd(Src(s).ang_in)];
%             Src(s).ID = 'I';
%             Src(s).AOI = 90;
%             Src(s).elem = pointLocation(mesh.tri,Src(s).Loc(1),Src(s).Loc(2));
%              Src(s).node = nearestNeighbor(mesh.tri,Src(s).Loc(:,1),Src(s).Loc(:,2));
%         elseif Src(s).Loc(1,2)==Domain(3)
%             %         normal_in = '+z'
%             Src(s).ang_in =temp; %180-Src(s).AOI ;
%             Src(s).Loc = Src(s).Loc + [mfp*sind(Src(s).ang_in),mfp*cosd(Src(s).ang_in)];
%             Src(s).ID = 'I';
%             Src(s).AOI = 90;
%             Src(s).elem = pointLocation(mesh.tri,Src(s).Loc(1),Src(s).Loc(2));
%             Src(s).node = nearestNeighbor(mesh.tri,Src(s).Loc(:,1),Src(s).Loc(:,2));
%         elseif Src(s).Loc(1,2)==Domain(4)
%             Src(s).ang_in=360-Src(s).AOI;
%             Src(s).Loc = Src(s).Loc + [mfp*sind(Src(s).ang_in),mfp*cosd(Src(s).ang_in)];
%             Src(s).ID = 'I';
%             Src(s).AOI = 90;
%             Src(s).elem = pointLocation(mesh.tri,Src(s).Loc(1),Src(s).Loc(2));
%              Src(s).node = nearestNeighbor(mesh.tri,Src(s).Loc(:,1),Src(s).Loc(:,2));
%         end
%         
%         get_one_mfp_source(mesh,Src(s),mfp)
%         
%     end
% end
%%

solver_options.pfuncflag = pfuncflag; 
solver_options.evalFluence = evalFluence; 
solver_options.probType = probType; 

if ~fcs_flag
    %basic spn
    [Jd]=func(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,solver_options);
else
    % fcs/ de variants
    [phi_u,Jd]=func(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,solver_options);
    varargout{1}=phi_u;
end
if (norm_flag==1)&&(simulation=='F')
    for s=1:NofSources
        if snr==100
            Jd_n=Jd{s}(:,SimType+1)./Jd{s}(:,SimType);
            if isreal(Jd_n)
                Md{s}=log(Jd_n);
            else
                Md.real{s}=real(log(Jd_n));
                Md.imag{s}=imag(log(Jd_n));
            end
        else
            Jd_n=awgn(Jd{s}(:,SimType+1)./Jd{s}(:,SimType),snr,'measured');
            if isreal(Jd_n)
                Md{s}=log(Jd_n);
            else
                Md.real{s}=real(log(Jd_n));
                Md.imag{s}=imag(log(Jd_n));
            end
        end
    end
else
    for s=1:NofSources
        if snr==100
            Jd_n=Jd{s}(:,SimType+1);
            if isreal(Jd_n)
                Md{s}=log(Jd_n);
            else
                Md.real{s}=real(log(Jd_n));
                Md.imag{s}=imag(log(Jd_n));
            end
        else
            Jd_n=awgn(Jd{s}(:,SimType+1),snr,'measured');
            if isreal(Jd_n)
                Md{s}=log(Jd_n);
            else
                Md.real{s}=real(log(Jd_n));
                Md.imag{s}=imag(log(Jd_n));
            end

%             if isreal(Jd_n)
%                 Md{s}=(Jd_n);
%             else
%                 Md.real{s}=real((Jd_n));
%                 Md.imag{s}=imag((Jd_n));
%             end

        end
    end
end
