% the idea is to compute a partial Jaccard index which could possibly disambiguate
% the role of a network C in the interaction between two other networks A
% and B. The influence of C can be uncorrelated, indirect, or synergistic
% uses the GRAMM toolbox https://github.com/piermorel/gramm
clear;clc;close all
n_nodes=50;
n_points=1000;
Jtot=zeros(1,3*n_points);
Jptot=zeros(1,3*n_points);


%% uncorrelated case

p1=0.1;
p2=0.1;
p3=0.1;

for m=1:n_points
    A=zeros(n_nodes);
    B=zeros(n_nodes);
    C=zeros(n_nodes);
    for i=1:n_nodes
        for j=i+1:n_nodes
            if rand<p1
                A(i,j)=1;
            end
            if rand<p2
                B(i,j)=1;
            end
            if rand<p3
                C(i,j)=1;
            end
        end
    end
    A=A+A';
    B=B+B';
    C=C+C';
    %let's compute Jaccard between A and B
    inda=find(A);
    indb=find(B);
    J(m)=length(intersect(inda,indb))/length(union(inda,indb));
    %partial Jaccard
    indc=find(C);
    inda=setdiff(inda,indc);indb=setdiff(indb,indc);
    Jp(m)=length(intersect(inda,indb))/length(union(inda,indb));
end
Jtot(1:n_points)=J;Jptot(1:n_points)=Jp;

%% indirect case

p1=0.1;
p2=0.1;
p3=0.1;

for m=1:n_points
    A=zeros(n_nodes);
    B=zeros(n_nodes);
    C=zeros(n_nodes);
    for i=1:n_nodes
        for j=i+1:n_nodes
            if rand<p1
                A(i,j)=1;
            end
            if rand<p2
                B(i,j)=1;
            end
            if rand<p3
                C(i,j)=1;
                A(i,j)=1;
                B(i,j)=1;
            end
        end
    end
    A=A+A';
    B=B+B';
    C=C+C';
    %Jaccard between A and B
    inda=find(A);
    indb=find(B);
    J(m)=length(intersect(inda,indb))/length(union(inda,indb));
    %partial Jaccard
    indc=find(C);
    inda=setdiff(inda,indc);indb=setdiff(indb,indc);
    Jp(m)=length(intersect(inda,indb))/length(union(inda,indb));
end
Jtot(n_points+1:2*n_points)=J;Jptot(n_points+1:2*n_points)=Jp;

%% suppression case

p1=0.1;
p2=0.1;
p3=0.1;
p4=1;
for m=1:n_points
    A=zeros(n_nodes);
    B=zeros(n_nodes);
    C=zeros(n_nodes);
    for i=1:n_nodes
        for j=i+1:n_nodes
            if rand<p1
                A(i,j)=1;
            end
            if rand<p2
                B(i,j)=1;
            end
            if rand<p3
                C(i,j)=1;
            end
            if A(i,j)*C(i,j)==0 && A(i,j)+C(i,j)>0
                if rand<p4
                    B(i,j)=1;
                end
            end
        end
    end
    A=A+A';
    B=B+B';
    C=C+C';
    %Jaccard between A e B
    inda=find(A);
    indb=find(B);
    J(m)=length(intersect(inda,indb))/length(union(inda,indb));
    %partial Jaccard
    indc=find(C);
    inda=setdiff(inda,indc);indb=setdiff(indb,indc);
    Jp(m)=length(intersect(inda,indb))/length(union(inda,indb));
end
Jtot(2*n_points+1:3*n_points)=J;Jptot(2*n_points+1:3*n_points)=Jp;

%% plot using Gramm

for i=1:n_points;lab{i}='1';end
for i=n_points+1:2*n_points;lab{i}='2';end
for i=2*n_points+1:3*n_points;lab{i}='3';end
clear g
g=gramm('x',Jtot,'y',Jptot,'color',lab);
g.geom_point();
%g.stat_cornerhist('edges',-.3:0.01:.5,'location',.5,'aspect',0.5,'fill','all');
g.set_text_options('base_size',14);
g.geom_abline();g.set_names('x','NJ(A,B)','y','NJ_{p}(A,B|C)','color','Algorithm');
g.set_text_options('base_size',14,'interpreter','tex');
g.set_color_options('map',[0 0 0.8; .5 .9 .2; 1 0 1]);
%g.axe_property('XLim',[.0 .5],'YLim',[.0 .5]);
g.draw();
%axis square
%set(gca,'LooseInset',get(gca,'TightInset'));
%set(gca,'Color','none');xlim([-0.01 1]);ylim([-0.01 1]);
