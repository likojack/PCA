clear all, close all

xin=[0.00 9.2830e-05 339.9 57.9
2.50 7.5820e-05 325.7 49.1
5.00 6.2920e-05 321.3 45.9
10.00 5.2090e-05 314.8 41.7
15.00 4.4550e-05 310.3 38.7
20.00 3.9540e-05 305.0 37.0
30.00 3.2570e-05 303.9 34.7
40.00 2.5670e-05 303.0 32.3
50.00 2.2520e-05 303.6 32.4
60.00 1.9820e-05 299.8 30.8
70.00 1.3890e-05 292.5 31.0
80.00 1.2570e-05 297.0 25.6
90.00 0.5030e-05 299.3 11.3];
%% col1: C,col2: AE, col3: AF, col4: AG
id0=find(xin(:,1)<16);
id1=find(xin(:,1)>16);
xyz=ID2XYZ(deg2rad(xin(:,4)),deg2rad(xin(:,3)));
xyz=bsxfun(@times,xyz,xin(:,2));

[coeff,score]=princomp(xyz(id1,:));
hat=[min(score(:,1));max(score(:,1))]*coeff(:,1)';
hat=bsxfun(@plus,hat,mean(xyz(id1,:)));


%% plot zplot with all points
h1=figure
subplot(2,2,1)
plot(xyz(:,2),xyz(:,1),'ok-','markerfacecolor','k')
hold on
plot(xyz(:,1),-xyz(:,3),'sk-','markerfacecolor','w')

%set(gca,'fontsize',12)
set(gca,'dataaspectratio',[1 1 1])
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
set(gca,'xlim',xlim,'ylim',ylim)
set(gca,'tickdir','out','xminortick','on','yminortick','on')
plot(xlim,[0 0],'--k')
plot([0 0],ylim,'--k')
%set(gca,'position',[0.1300    0.5838    0.3347    0.3412])
set(gca,'ticklength',get(gca,'ticklength')*1.5)

subplot(2,2,2)
plot(xyz(id1,2),xyz(id1,1),'ok','markerfacecolor','k')
hold on
plot(xyz(id1,1),-xyz(id1,3),'sk','markerfacecolor','w')
plot(hat(:,2),hat(:,1),'r','linewidth',1.5)
plot(hat(:,1),-hat(:,3),'r','linewidth',1.5)

%set(gca,'fontsize',10)
set(gca,'dataaspectratio',[1 1 1])
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
set(gca,'xlim',xlim,'ylim',ylim)
set(gca,'tickdir','out','xminortick','on','yminortick','on')
plot(xlim,[0 0],'--k')
plot([0 0],ylim,'--k')
set(gcf,'position',[680         129        1118         849])
set(gca,'ticklength',get(gca,'ticklength')*1.5)
xtick=get(gca,'xtick')