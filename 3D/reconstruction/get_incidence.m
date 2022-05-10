function Incid_mat = get_incidence(DT)
NofElem = size(DT,1);
Ne=neighbors(DT);
elem_nums = repmat((1:NofElem)',[3,1]);
graph_edge(:,1)=elem_nums(:);
graph_edge(1:NofElem,2)=Ne(:,1);
graph_edge(NofElem+1:2*NofElem,2)=Ne(:,2);
graph_edge(2*NofElem+1:3*NofElem,2)=Ne(:,3);

graph_edge_sort=unique(sort(graph_edge,2),'rows');
graph_final=[];
temp =find(isnan(graph_edge_sort(:,2)));
last_ind=0;
for ii=1:length(temp)
    graph_final=[graph_final;graph_edge_sort(last_ind+1:temp(ii)-1,:)];
    last_ind=temp(ii);
end
Incid_mat=zeros(length(graph_final),NofElem);
for ii=1:length(graph_final)
    Incid_mat(ii,graph_final(ii,1))=-1;
    Incid_mat(ii,graph_final(ii,2))=1;
end