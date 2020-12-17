clear;clc;load('copenhagen_new.mat')
%%%% copenhagen

N_nodes=851;
% 1 fb,CS,BT; 2 BT,CS,fb; 3 BT,fb,CS
%Distance proxy - Received Signal Strength (dBm)

%load Bnew
myB=BTW1_far;
%myB=BTW7_far;%col=c(7,:);
%CS=CS(setdiff(1:N_nodes,discon),setdiff(1:N_nodes,discon));
%fb=fb(setdiff(1:N_nodes,discon),setdiff(1:N_nodes,discon));
A1=myB;
A2=fb;
A3=CS; 
N_nodes=length(A1);
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