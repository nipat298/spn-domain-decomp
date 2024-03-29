function partmat=OptimalSchwarz(partmat,alpha,N)
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
% max_tol=max([partmat(:).tol]);

iter=0;
% while iter>0
 while max_tol>0.005 || iter<=2
    parfor i=1:NoDomains
        len_node=length(partmat(i).nodes);
        % The RHS_temp variable is the RHS_extra variable from noMETIS_4sub code
        RHS_temp=zeros(((N+1)/2)*len_node,1);
        for j=1:NoDomains
            if j==i
                continue
            else
               int_node=partmat(i).int_local{j};
               len_sp3= int_node(:,1) + len_node;
               partmat(i).robin{j}(int_node)=g_ini{i,j}(1:length(int_node));
               partmat(i).robin{j}(len_sp3)=g_ini{i,j}((1+length(int_node)):end);
                
%                partmat(i).robin{j}(int_node+len_node)=g_ini{i,j}((1+length(int_node)):end);
            
               RHS_temp=RHS_temp+(partmat(i).RHS_extra{j}*partmat(i).robin{j});
            end
        end
        partmat(i).RHS_final=partmat(i).source+RHS_temp;
        partmat(i).soln=gmres(partmat(i).Gsub,partmat(i).RHS_final,[],tol,maxit);
%     partmat(i).soln=partmat(i).Gsub\partmat(i).RHS_final;
        partmat(i).tol=abs(partmat(i).soln-partmat(i).uprev)./abs(partmat(i).soln);
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
                u_int{i,j}((1:length(int_node)),1)=partmat(i).soln(int_node);
                u_int{i,j}(((1+length(int_node)):(2*length(int_node))),1)=partmat(i).soln(len_sp3);
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
        partmat(i).uprev(:,:)=partmat(i).soln(:,:);
        partmat(i).maxtol=max(partmat(i).tol);
    end
    
    g_ini(:,:)=g_next(:,:);
  max_tol=max([partmat.maxtol]);
    iter=iter+1;
    disp(iter);
%      iter=iter-1;


toc
 end