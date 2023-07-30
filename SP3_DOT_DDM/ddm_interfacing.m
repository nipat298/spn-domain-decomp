
function partmat= ddm_interfacing(mesh,NoDomains,varargin)
% This function performs subdomain partitioning along the axis provided by the user.
% It also classifies nodes into respective subdomains and lists elements across original boundary and interface.

%   Inputs:
%   mesh: structure type variable which stores meshing information
%   NoDomains: Number of subdomains
%   varargin: Variable arguments for the current code include int_axes and Domain

%   Return variables:
%   partmat: Structure type variable which contains information for all subdomains

int_axes=varargin{1};
Domain=varargin{2};

ElemConnMatrix = [mesh.tri(:,1),mesh.tri(:,2),mesh.tri(:,3),mesh.tri(:,4)];
NofNodes=size(mesh.nodes,1);


t=ElemConnMatrix(:,:);
Nfaces=size(mesh.tri,2);
Nelem=size(mesh.tri,1);

elem=(1:Nelem)';
elem=repmat(elem,1,Nfaces);

% This switch case helps helps in partitioning along the choden axis
switch(int_axes)
    case 'x'
        dmin=Domain(1);
        dmax=Domain(2);
        coord=1;
    case 'y'
        dmin=Domain(3);
        dmax=Domain(4);
        coord=2;
    case 'z'
        dmin=Domain(5);
        dmax=Domain(6);
        coord=3;
end
sd_len=(dmax-dmin)/NoDomains;
tol = 1e-6; % Tolerance to be used with 'find' function

% This section gives interface planes and interface nodes
% Here interface is also a structure type variable with fields giving interface planes and nodes along each interface plane
interface_nodes=[]; % This variable stores nodes across all the interfaces combined
for kk=1:(NoDomains-1)
    interface(kk).planes= dmin + kk*sd_len;
    disp(interface(kk).planes);
    nod=find((abs(mesh.nodes(:,coord)-interface(kk).planes))>=0 & abs((mesh.nodes(:,coord)-interface(kk).planes))<tol);
    interface(kk).nodes=nod;
    interface_nodes=[interface_nodes interface(kk).nodes];
end

for kk=1:NoDomains
    Subdomain(kk).nodes=[];
end

% This section gives nodes in each subdomain
for cc=1:NofNodes
    for kk=1:NoDomains
        if kk==1
            if mesh.nodes(cc,coord)>=dmin && mesh.nodes(cc,coord)<interface(kk).planes
                Subdomain(kk).nodes=[Subdomain(kk).nodes;cc];
            end
        elseif kk>1 && kk<NoDomains
            if mesh.nodes(cc,coord)>interface(kk-1).planes && mesh.nodes(cc,coord)<interface(kk).planes
                Subdomain(kk).nodes=[Subdomain(kk).nodes;cc];
            end
        elseif kk==NoDomains
            if mesh.nodes(cc,coord)>interface(kk-1).planes && mesh.nodes(cc,coord)<=dmax
                Subdomain(kk).nodes=[Subdomain(kk).nodes;cc];
            end
        end
    end
end

% Subdomain is a structure type variable with the fields representing nodes, global interface numbering for each subdomain

% This section gives global numbering for interface nodes of each subdomain
for kk=1:NoDomains
    if kk==1
        Subdomain(kk).nodes=unique([Subdomain(kk).nodes; interface(kk).nodes]);
        Subdomain(kk).interface_g{kk+1}=interface(kk).nodes;
    elseif kk>1 && kk<NoDomains
        Subdomain(kk).nodes=unique([Subdomain(kk).nodes;interface(kk-1).nodes;interface(kk).nodes]);
        Subdomain(kk).interface_g{kk-1}=interface(kk-1).nodes;
        Subdomain(kk).interface_g{kk+1}=interface(kk).nodes;
    elseif kk==NoDomains
        Subdomain(kk).nodes=unique([Subdomain(kk).nodes;interface(kk-1).nodes]);
        Subdomain(kk).interface_g{kk-1}=interface(kk-1).nodes;
    end
end

% This section transfers the data stored in the structure variable 'Subdomain' to the structure variable 'partmat'
for kk=1:NoDomains
    partmat(kk).nodes=Subdomain(kk).nodes;
    partmat(kk).interface_g=Subdomain(kk).interface_g;
end


for i=1:NoDomains
    for j=1:NoDomains
        if j==(i-1) || j==(i+1)
            continue
        else 
            partmat(i).interface_g{j}=[0];
        end
    end
end

% This section lists out the total number of interface nodes for each subdomain and stores them as a new field in partmat variable
for i=1:NoDomains
    if i==1
        partmat(i).inter_tot=interface_nodes(:,i);
    elseif i==NoDomains
        partmat(i).inter_tot=interface_nodes(:,(i-1));
    elseif i>1 && i<NoDomains
        partmat(i).inter_tot=unique([interface_nodes(:,(i-1));interface_nodes(:,i)]);
    end
end

% This section lists out the global node numbers of nodes in each subdomain
for i=1:NoDomains
    NC_sub=length(partmat(i).nodes);
    partmat(i).SD_nodes=[(1:NC_sub)' partmat(i).nodes];
end

% This section lists out the local node numbers for each subdomain
for i=1:NoDomains
    for j=1:NoDomains
            ismem=ismember(partmat(i).nodes,partmat(i).interface_g{j});
            partmat(i).int_local{j}=find(ismem);
    end
end

% This function forms Element Connectivity Matrix for each subdomain
partmat=ECM_3D(partmat,ElemConnMatrix,NoDomains);
 
%Classifying original boundary triangles and interface triangles for each subdomain

partmat=boundtri_bdelem_3D(partmat,NoDomains); % To obtain boundary triangles for each subdomain
for i=1:NoDomains
    tt1=ismember(partmat(i).boundtri,partmat(i).inter_tot);
    tt2=find(tt1(:,1)==1 & tt1(:,2)==1 & tt1(:,3)==1);
    
    % Following fields give triangles and elements along interface for each subdomain
    partmat(i).art_bound=partmat(i).boundtri(tt2,:);
    partmat(i).artb_elem=partmat(i).belem(tt2,:);
    
    % Following fields give triangles and elements along original boundary for each subdomain
    partmat(i).boundtri(tt2,:)=[];
    partmat(i).belem(tt2,:)=[];
    
end
