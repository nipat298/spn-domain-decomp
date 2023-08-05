function partmat=OptimalSchwarz_ex(partmat,lambda,N)
% This function implements the domain decomposition algorithm for excitation equation

%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partmat: Structure type variable containing all the pre processing information for each subdomain
% lambda: Interface parameter (transmission coefficient)
% N: order of the SPN approximation (3 for SP3 case)

%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%
% partmat: Structure type variable with the final solution across all subdomains included

% For better understanding of this subroutine, refer the text file "Detailed code layout" shared with the SP3_DOT_DDM folder
tic
tol = 1e-6;
maxit = 1000;
NoDomains=size(partmat,2);

g_ini=cell(NoDomains,NoDomains); % Initial assumption for interface nodes

for i=1:NoDomains
    len_sub=((N+1)/2)*length(partmat(i).nodes);
    for j=1:NoDomains
        if j==i
            continue
        else
            len_int=length(partmat(i).int_local{j});
             % g_int is the initial value we are assuming for the interface nodes. 
            % Here we have assumed an initial "zero" value for all interface nodes
            g_int=zeros(((N+1)/2)*len_int,1);
            % SP3 is a coupled system of equations with two composite moments.
            % Therefore if an interface has N nodes, then for SP3 case, we define a variable with the size 2N for that interface.
            % 1:N will correspond to 1st composite moment 
            % (N+1):2N will correspond to 2nd composite moment
            g_ini{i,j}=g_int;
            partmat(i).robin{j}=zeros(len_sub,1);
        end
    end
    
    partmat(i).uprev=zeros(len_sub,1);
    partmat(i).tol=1;
end

max_tol=1;
 %max_tol is the convergence criteria (defined in section 3.1 and equation (3.2) of ref 3)

iter=0;
while max_tol>0.01 || iter<=2
    parfor i=1:NoDomains
        len_node=length(partmat(i).nodes);
        RHS_temp=zeros(((N+1)/2)*len_node,1);
        for j=1:NoDomains
            if j==i
                continue
            else
               int_node=partmat(i).int_local{j};
               len_sp3= int_node(:,1) + len_node;
               partmat(i).robin{j}(int_node)=g_ini{i,j}(1:length(int_node));
               partmat(i).robin{j}(len_sp3)=g_ini{i,j}((1+length(int_node)):end);
                % This calculation is done for each interface of a subdomain
               % RHS_temp is the extra term that is added to RHS because of DDM algorithm
               RHS_temp=RHS_temp+(partmat(i).RHS_extra{j}*partmat(i).robin{j});
            end
        end
        partmat(i).RHS_final=partmat(i).source+RHS_temp;
        partmat(i).soln=gmres(partmat(i).Gsub,partmat(i).RHS_final,[],tol,maxit);% Using GMRES solver to obtain solution
        partmat(i).tol=abs(partmat(i).soln-partmat(i).uprev)./abs(partmat(i).soln); % Finding relative error using current and previous iterate
    end
    
   % This section extracts solution value at interface nodes for each subdomain
   % and stores them in a cell form of structure 'u_int'. This section runs for all
   % interfaces of a given subdomain
    u_int=cell(NoDomains,NoDomains);
    for i=1:NoDomains
        len_node=length(partmat(i).nodes);
        for j=1:NoDomains
            if j==i
                continue
            else
                int_node=partmat(i).int_local{j};
                len_sp3= int_node(:,1) + len_node;
                u_int{i,j}((1:length(int_node)),1)=partmat(i).soln(int_node);% For 1st composite moment
                u_int{i,j}(((1+length(int_node)):(2*length(int_node))),1)=partmat(i).soln(len_sp3);% For 2nd composite moment
            end
        end
    end
    
    % The following section handles the update condition for DDM algorithm.
    % This update condition is applied after each iteration
    g_next=cell(NoDomains,NoDomains);% Just like g_ini, we initialize a cell structure for g_next
    for i=1:NoDomains
        for j=1:NoDomains
            if j==i
                continue
            else 
                g_next{i,j}=2*lambda*u_int{j,i}-g_ini{j,i};% This is the update condition for DDM algorithm
            end
        end
         %Updating the value of 'uprev' vector with the current iterate solution vector
        partmat(i).uprev(:,:)=partmat(i).soln(:,:);
        partmat(i).maxtol=max(partmat(i).tol);
    end
    
    g_ini(:,:)=g_next(:,:);%Updating g_ini for next iteration
  max_tol=max([partmat.maxtol]);
    iter=iter+1;%Counter for iteration number
    disp(iter);

toc
 end