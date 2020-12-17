clear;clc
%importNYC;
load NYC
%%%%  1 retweet 2 mentions 3 replies
% 1 {1,2,3} 2 {1,3,2} 3 {2,3,1}
N_nodes=102439;
i1=find(MP(:,1)==2);
i2=find(MP(:,1)==3);
L1=MP(i1,2:4);
L2=MP(i2,2:4);
A1=sparse(N_nodes); % for small matrices you can use zeros() instead of sparse()
for i=1:length(L1);A1(L1(i,1),L1(i,2))=L1(i,3);end
A2=sparse(N_nodes);
for i=1:length(L2);A2(L2(i,1),L2(i,2))=L2(i,3);end
A1(N_nodes,N_nodes)=0;A2(N_nodes,N_nodes)=0;
AA2=A2;AA1=A1;
A1=reshape(A1,N_nodes*N_nodes,1);
A2=reshape(A2,N_nodes*N_nodes,1);
i1=find(A1);
i2=find(A2);
J=length(intersect(i1,i2))/length(union(i1,i2));
%partial Jaccard

% random shuffling
parfor h=1:500
    
    
    
    % complete rewiring
    %     jj=randperm(N_nodes);
    %     A1=reshape(AA1(jj,jj),N_nodes*N_nodes,1);
    %     jj=randperm(N_nodes);
    %     A2=reshape(AA2(jj,jj),N_nodes*N_nodes,1);
    jj=randperm(N_nodes);
    A1=reshape(AA1(jj,jj),N_nodes*N_nodes,1);
        jj=randperm(N_nodes);
    A2=reshape(AA2(jj,jj),N_nodes*N_nodes,1);
    ind1=find(A1);
    ind2=find(A2);
    Js1(h)=length(intersect(ind1,ind2))/length(union(ind1,ind2));   
end
fprintf('%.7f %.7f %.7f\n',J,mean(Js1),std(Js1));