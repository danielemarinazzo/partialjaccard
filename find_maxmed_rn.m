function maxmed=find_maxmed_rn(A1,A2,A3)
AA3=A3;
N_nodes=length(A1);
A1=reshape(A1,N_nodes*N_nodes,1);
A2=reshape(A2,N_nodes*N_nodes,1);
A3=reshape(A3,N_nodes*N_nodes,1);
[r3]=find(A3);
% fA1=find(A1);
% lf=length(fA1);
% rf=randperm(lf);
% fA1r=fA1(rf);
% if length(r3)<lf
%     A1(fA1r(1:length(r3)))=0;
% end
% fA2=find(A2);
% lf=length(fA2);
% rf=randperm(lf);
% fA2r=fA2(rf);
% if length(r3)<lf
%     A2(fA2r(1:length(r3)))=0;
% end
A1(r3)=1;
A2(r3)=1;
i1=find(A1);
i2=find(A2);
%partial Jaccard
i3=find(A3);
i13=setdiff(i1,i3); %links which are in A1 but not in A3
i23=setdiff(i2,i3); %links which are in A2 but not in A3
 
i2_3=setdiff(i2,i23);% links which are both in A2 and A3
i231=setdiff(i23,i13);%links which are in A2 but not in A1 neither in A3
n=length(i23)-length(i231);%number of links which are both in A2 and A1, but not in A3
%in3=setdiff(find(ones(N_nodes)),i3);%links which are not in A3
i111=intersect(i1,i2);i111=intersect(i111,i3);%links common to A1, A2, and A3



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
    % remove synergy
%     i10syn=setdiff(i1,i3);i20syn=setdiff(i2,i3);
%     %i11=setdiff(i1,i10);
%     i21syn=setdiff(i2,i20syn);
%     i200syn=setdiff(i20syn,i10syn);
%     n=length(i20syn)-length(i200syn);
%     iisyn=setdiff(find(ones(N_nodes)),i3);
    %g=randperm(length(in3));
    rl=ceil(rand(n,1)*N_nodes^2); % check for intersect(rl,i3)
    i231syn=union(i231,rl); % adding random links instead of those for which A1=A2=1 and A3=0
    ind2syn=union(i231syn,i2_3); %adding links with A2=A3=1
    Js(h)=length(intersect(i1,ind2syn))/length(union(i1,ind2syn));
    %partial Jaccard
    Jps(h)=length(intersect(i13,i231syn))/length(union(i13,i231syn));
    
    % remove synergy on randomized
    i10sr=setdiff(i1,ind3);%links in A1 but not in the new A3 
    i20sr=setdiff(i2,ind3);%links in A2 but not in the new A3
    i11=setdiff(i1,i10sr);%links in A1 and A3
    i21=setdiff(i2,i20sr);%links in A2 and A3
    i200=setdiff(i20sr,i10sr);%links in A2 but not in A1 neither in A3
    nsr=length(i20sr)-length(i200);%number of links in A1 and A2 but not in A3
    %ii=setdiff(find(ones(N_nodes)),ind3);%
    %g=randperm(length(ii));
    rl=ceil(rand(nsr,1)*N_nodes^2); %check intersect(rl,ind3);
    i2000=union(i200,rl);% adding random links instead of those for which A1=A2=1 and A3=0
    ind2=union(i2000,i21);%adding links in A1 and A3
    Jsr(h)=length(intersect(i1,ind2))/length(union(i1,ind2));
    %     %partial Jaccard
    Jpsr(h)=length(intersect(i10sr,i2000))/length(union(i10sr,i2000));   
end

%%
med_index=(Jps-Js)-(Jpsr-Jsr);
%med_index(med_index<0)=0;
maxmed=mean(med_index);