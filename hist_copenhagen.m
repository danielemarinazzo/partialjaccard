load copenhagen_new
c=autumn(7);
figure;
%[N,edges] = histcounts(log10(-Bnonzero),20);
[N,edges] = histcounts(Bnonzero,40);
B1=BT1_far(BT1_far~=0);
hold on;area([(min(B1)) (max(B1))],[max(N) max(N)],'FaceAlpha',1,'FaceColor',c(7,:))
B2=BT2_far(BT2_far~=0);
hold on;area([(min(B2)) (max(B2))],[max(N) max(N)],'FaceAlpha',1,'FaceColor',c(6,:))
B3=BT3_far(BT3_far~=0);
hold on;area([(min(B3)) (max(B3))],[max(N) max(N)],'FaceAlpha',1,'FaceColor',c(5,:))
B4=BT4_far(BT4_far~=0);
hold on;area([(min(B4)) (max(B4))],[max(N) max(N)],'FaceAlpha',1,'FaceColor',c(4,:))
B5=BT5_far(BT5_far~=0);
hold on;area([(min(B5)) (max(B5))],[max(N) max(N)],'FaceAlpha',1,'FaceColor',c(3,:))
B6=BT6_far(BT6_far~=0);
hold on;area([(min(B6)) (max(B6))],[max(N) max(N)],'FaceAlpha',1,'FaceColor',c(2,:))
B7=BT7_far(BT7_far~=0);
hold on;area([(min(B7)) (max(B7))],[max(N) max(N)],'FaceAlpha',1,'FaceColor',c(1,:))
%hold on;histogram(-Bnonzero,10.^edges,'FaceColor','k','EdgeColor','none','FaceAlpha',.6);set(gca,'xscale','log')
hold on;histogram(Bnonzero,edges,'FaceColor','k','EdgeColor','none','FaceAlpha',.6);%set(gca,'xscale','log')
set(gca,'YTick',[])
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'));
%set(gca,'Color','none');