
%%% This function gives the ElemConnMatrix for each subdomain.

%%%% Inputs:
% partmat: Struct containing nodes info of each subdomain. It's a struct for the subdomains and their different parameters
% EkemConnMatrix: The global ElemConnMatrix for the entire domain
% N: number of subdomains

%%%% Output:
% partmat: Struct with an added field for ECM of each subdomain

function partmat=ECM_3D(partmat,ElemConnMatrix,N)
for i=1:N
    ind=ismember(ElemConnMatrix,partmat(i).nodes);
    elem_index=ind(:,1)&ind(:,2)&ind(:,3)&ind(:,4);
    elem=find(elem_index);
    ECM_SD=ElemConnMatrix(elem,:);
    partmat(i).ECM=[ECM_SD elem];
end
end

