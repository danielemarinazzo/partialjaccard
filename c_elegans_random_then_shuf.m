clear;clc;load('celegans.mat')
 
%%%% c-elegans 1 electric 2 chem monadic 3 chem polyadic
% 1 {1,2,3} 2 {1,3,2} 3 {3,2,1}
N_nodes=279;
col=[230 159 0]/255;
i1=find(MP(:,1)==3);
i2=find(MP(:,1)==2);
i3=find(MP(:,1)==1);
L1=MP(i1,2:4);
L2=MP(i2,2:4);
L3=MP(i3,2:4);
A1=zeros(N_nodes); % for small matrices you can use zeros() instead of sparse()
for i=1:length(L1);A1(L1(i,1),L1(i,2))=L1(i,3);end
A2=zeros(N_nodes);
for i=1:length(L2);A2(L2(i,1),L2(i,2))=L2(i,3);end
A3=zeros(N_nodes);
for i=1:length(L3);A3(L3(i,1),L3(i,2))=L3(i,3);end
AA3=A3;AA2=A2;AA1=A1;
A1=reshape(A1,N_nodes*N_nodes,1);
A2=reshape(A2,N_nodes*N_nodes,1);
A3=reshape(A3,N_nodes*N_nodes,1);
i1=find(A1);
i2=find(A2);
J=length(intersect(i1,i2))/length(union(i1,i2));
%partial Jaccard
i3=find(A3);
i13=setdiff(i1,i3); %links which are in A1 but not in A3
i23=setdiff(i2,i3); %links which are in A2 but not in A3
Jp=length(intersect(i13,i23))/length(union(i13,i23));
 
i2_3=setdiff(i2,i23);% links which are both in A2 and A3
i231=setdiff(i23,i13);%links which are in A2 but not in A1 neither in A3
n=length(i23)-length(i231);%number of links which are both in A2 and A1, but not in A3
in3=setdiff(find(ones(N_nodes)),i3);%links which are not in A3
i111=intersect(i1,i2);i111=intersect(i111,i3);%links common to A1, A2, and A3
iu=setdiff(i1,i111);%links which are in A1 but not in the intersection of A2 and A3
n111=length(i111);%number of links common to A1, A2, and A3
 
 
% random shuffling
parfor h=1:500
    
    
    
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
figure
h1=histogram(Jps1-J,'FaceColor','g','FaceAlpha',.5);
hold on
h2=histogram(Jps-Js,'FaceColor',[0 51 255]/255,'FaceAlpha',.5);
h3=histogram(Jpr-Jr,'FaceColor',[0.5430 0 0],'FaceAlpha',.5);
h4=histogram(Jpsr-Jsr,'FaceColor',[102 204 255]/255,'FaceAlpha',.5);
h5=histogram(Jprr-Jrr,'FaceColor',[0.9792 0.5000 0.4453],'FaceAlpha',.5);
val=[max(h1.Values);max(h2.Values);max(h3.Values)];%;max(h4.Values);max(h5.Values)];
mval=max(val);
close
% figure;
% %axes('Color','none','YColor','none');
% line([Jp-J Jp-J],[0 mval],'Color','k','LineWidth',2);hold on;
% h1=histogram(Jps1-J,'FaceColor',[0 0 0.8],'FaceAlpha',.5,'EdgeColor','none');
% h2=histogram(Jps-Js,'FaceColor',[.5 .9 .2],'FaceAlpha',1,'EdgeColor','none');
% h3=histogram(Jpr-Jr,'FaceColor',[1 0 1],'FaceAlpha',1,'EdgeColor','none');
% h4=histogram(Jpsr-Jsr,'FaceColor',[.9 1 .2],'FaceAlpha',.5,'EdgeColor','none');
% h5=histogram(Jprr-Jrr,'FaceColor',[1 .7 .8],'FaceAlpha',.5,'EdgeColor','none');box on
% val=[max(h1.Values);max(h2.Values);max(h3.Values)];%;max(h4.Values);max(h5.Values)];
% mval=max(val);
% % legend('no rewiring','rewiring of C', 'remove suppression', 'remove mediation',...
% %     'remove suppression after rewiring', 'remove mediation after rewiring','location','westoutside');
% legend('\Delta_{0}','R', 'S', 'M','S,R', 'M,R');%,'location','westoutside');
% legend boxoff
% %title('C Elegans: Electric, Chemical Monadic, Chemical Polyadic')
% xlabel('\Delta')
% set(gca, 'FontSize',20,'YTick',[])
% set(gca,'LooseInset',get(gca,'TightInset'));
% set(gcf,'Position',[10 10 600 500])
%% 
Delta=Jp-J;
Delta_M=Jps-Js;
Delta_RM=Jpsr-Jsr;
Delta_S=Jpr-Jr;
Delta_RS=Jprr-Jrr;
Delta_R=Jps1-J;
for i=1:20
%MM(i)=find_maxmed(N_nodes,round(length(i3)/2));
%MS(i)=find_maxsyn(N_nodes,round(length(i3)/2));
MM(i)=find_maxmed_rn(AA1,AA2,AA3);
MS(i)=find_maxsyn_rn(AA1,AA2,AA3);
end
maxmed=mean(MM);maxsyn=mean(MS);
syn_index=(Delta_M-mean(Delta_RM))./maxsyn;
med_index=-(Delta_S-mean(Delta_RS))./maxmed;
syn_sig=abs((Delta_M-mean(Delta_RM))./std(Delta_RM));
med_sig=abs((Delta_S-mean(Delta_RS))./std(Delta_RS));
%med_index=-(Delta_RM-Delta_M);
%syn_index=-(Delta_RS-Delta_S);
%med_sig=(Delta_RM-mean(Delta_M))/std((Delta_RM-mean(Delta_M)));
%syn_sig=(Delta_RS-mean(Delta_S))/std((Delta_RS-mean(Delta_S)));
%med_index(med_index<0)=0;
%syn_index(syn_index>0)=0;syn_index=-syn_index;
h1=figure(5);hold on;scatter(syn_index,med_index,25,col,'filled');
xlabel('\Delta_{W,S}');ylabel('\Delta_{W,M}');set(gca,'FontSize',30)
axis square
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'Color','none');%xlim([-0.01 1]);ylim([-0.01 1]);
h1=figure(6);hold on;scatter(syn_sig,med_sig,25,col,'filled');
xlabel('\Delta_{W,S}');ylabel('\Delta_{W,M}');set(gca,'FontSize',30)
axis square
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'Color','none');%xlim([-0.01 1]);ylim([-0.01 1]);