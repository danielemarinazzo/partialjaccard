clear;clc;load('celegans.mat')
 
%%%% c-elegans 1 electric 2 chem monadic 3 chem polyadic
% 1 {3,2,1} 2 {1,3,2} 3 {1,2,3}
nmaps=3;c=colormap(lines(nmaps));
N_nodes=279;
col=c(2,:);
i1=find(MP(:,1)==1);
i2=find(MP(:,1)==2);
L1=MP(i1,2:4);
L2=MP(i2,2:4);
A1=zeros(N_nodes); % for small matrices you can use zeros() instead of sparse()
for i=1:length(L1);A1(L1(i,1),L1(i,2))=L1(i,3);end
A2=zeros(N_nodes);
for i=1:length(L2);A2(L2(i,1),L2(i,2))=L2(i,3);end
AA2=A2;AA1=A1;
A1=reshape(A1,N_nodes*N_nodes,1);
A2=reshape(A2,N_nodes*N_nodes,1);
i1=find(A1);
i2=find(A2);
J=length(intersect(i1,i2))/length(union(i1,i2));

AR1=zeros(size(A1));AR2=AR1; 
% random shuffling
parfor h=1:500
    jj=randperm(N_nodes^2);
    AR1=A1(jj);
    jj=randperm(N_nodes^2);
    AR2=A2(jj);
    ind1=find(AR1);
    ind2=find(AR2);
    Js1(h)=length(intersect(ind1,ind2))/length(union(ind1,ind2));
    
end
fprintf('%.7f %.7f %.7f\n',J,mean(Js1),std(Js1));