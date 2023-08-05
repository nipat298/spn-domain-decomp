function partmat=OptimalSchwarz_em(partmat,alpha,N)
% This function implements the domain decomposition algorithm for emission equation

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
            g_int=zeros(((N+1)/2)*len_int,1);
            g_ini{i,j}=g_int;
            partmat(i).robin{j}=zeros(len_sub,1);
        end
    end
    
    partmat(i).uprev=zeros(len_sub,1);
    partmat(i).tol=1;
end

max_tol=1;

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
               RHS_temp=RHS_temp+(partmat(i).RHS_extra{j}*partmat(i).robin{j});
            end
        end            
            
  
        %Replacing partmat(i).source by partmat(i).SM here
        partmat(i).RHS_final=partmat(i).SM+RHS_temp;
        partmat(i).solnFOT=gmres(partmat(i).Gsub,partmat(i).RHS_final,[],tol,maxit);
        partmat(i).tol=abs(partmat(i).solnFOT-partmat(i).uprev)./abs(partmat(i).solnFOT);
    end
    
   
    u_int=cell(NoDomains,NoDomains);
    for i=1:NoDomains
        len_node=length(partmat(i).nodes);
        for j=1:NoDomains
            if j==i
                continue
            else
                int_node=partmat(i).int_local{j};
                len_sp3= int_node(:,1) + len_node;
                u_int{i,j}((1:length(int_node)),1)=partmat(i).solnFOT(int_node);
                u_int{i,j}(((1+length(int_node)):(2*length(int_node))),1)=partmat(i).solnFOT(len_sp3);
            end
        end
    end
    
    g_next=cell(NoDomains,NoDomains);
    for i=1:NoDomains
        for j=1:NoDomains
            if j==i
                continue
            else 
                g_next{i,j}=2*alpha*u_int{j,i}-g_ini{j,i};
            end
        end
        partmat(i).uprev(:,:)=partmat(i).solnFOT(:,:);
        partmat(i).maxtol=max(partmat(i).tol);
    end
    
    g_ini(:,:)=g_next(:,:);
  max_tol=max([partmat.maxtol]);
    iter=iter+1;
    disp(iter);


toc
 end