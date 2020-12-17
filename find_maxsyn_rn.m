function maxsyn=find_maxsyn_rn(A1,A2,A3)
AA3=A3;
N_nodes=length(A1);
A1=reshape(A1,N_nodes*N_nodes,1);
A2=reshape(A2,N_nodes*N_nodes,1);
A3=reshape(A3,N_nodes*N_nodes,1);
P=A1.*A3;
S=A1+A3;
S(S>0)=1;
[r]=find(S+P==1);
% fA2=find(A2);
% lf=length(fA2);
% rf=randperm(lf);
% fA2r=fA2(rf);
% if length(r)<lf
%     A2(fA2r(1:length(r)))=0;
% end
A2(r)=1;
i1=find(A1);
i2=find(A2);


i3=find(A3);
i13=setdiff(i1,i3); %links which are in A1 but not in A3
i23=setdiff(i2,i3); %links which are in A2 but not in A3
 
i231=setdiff(i23,i13);%links which are in A2 but not in A1 neither in A3
%in3=setdiff(find(ones(N_nodes)),i3);%links which are not in A3
i111=intersect(i1,i2);i111=intersect(i111,i3);%links common to A1, A2, and A3
iu=setdiff(i1,i111);%links which are in A1 but not in the intersection of A2 and A3
n111=length(i111);%number of links common to A1, A2, and A3


% random shuffling
parfor h=1:20
    
        
    % complete rewiring
    %     jj=randperm(N_nodes);
    %     A1=reshape(AA1(jj,jj),N_nodes*N_nodes,1);
    %     jj=randperm(N_nodes);
    %     A2=reshape(AA2(jj,jj),N_nodes*N_nodes,1);
    jj=randperm(N_nodes);
    A3=reshape(AA3(jj,jj),N_nodes*N_nodes,1);
    %ind1=find(A1);
    %ind2=find(A2);
      ind1=i1;
      ind2=i2;
    %Js1(h)=length(intersect(ind1,ind2))/length(union(ind1,ind2));
    %partial Jaccard
    ind3=find(A3); %links in the new A3
    ind1=setdiff(ind1,ind3);ind2=setdiff(ind2,ind3); %subtraction of the A3 links from those of A1 and A2
    Jps1(h)=length(intersect(ind1,ind2))/length(union(ind1,ind2));
    % removing redundancy
%     i1r=intersect(i1,i2);i1r=intersect(i1r,i3);
%     iu=setdiff(i1,i1r);
%     nr=length(i1r);
    %g=randperm(N_nodes^2);
    rl=ceil(rand(n111,1)*N_nodes^2);
    inda1=union(iu,rl);%adding random links instead of those common to A1, A2 and A3
    indaa=setdiff(inda1,i3);
    Jr(h)=length(intersect(inda1,i2))/length(union(inda1,i2));
    %partial Jaccard
    Jpr(h)=length(intersect(indaa,i23))/length(union(indaa,i23));
    
    % removing redundancy from randomized
    i1r=intersect(i1,i2);i1r=intersect(i1r,ind3);
    iurr=setdiff(i1,i1r);
    nrr=length(i1r);
    %g=randperm(N_nodes^2);
    rl=ceil(rand(nrr,1)*N_nodes^2);
    inda1=union(iurr,rl);
    indaa=setdiff(inda1,ind3);
    i23n=setdiff(i2,ind3);
    Jrr(h)=length(intersect(inda1,i2))/length(union(inda1,i2));
    %partial Jaccard
    Jprr(h)=length(intersect(indaa,i23n))/length(union(indaa,i23n)); 
end


%%
syn_index=-((Jpr-Jr)-(Jprr-Jrr));
%syn_index(syn_index>0)=0;syn_index=-syn_index;
maxsyn=mean(syn_index);
