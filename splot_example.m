%Example of how to produce a stereonet

%% equal area projection
figure
t=linspace(0,2*pi)';
h=fill(sin(t),cos(t),'w')
axis square, hold on

%plot circles representing constant inclinations
I0=[15:15:75]';
for i=1:numel(I0)
    D1=[0:1:360]';
    I1=ones(size(D1)).*I0(i);
    XYZtick=ID2XYZ(deg2rad(I1),deg2rad(D1));
    Stick=pmag_splot(XYZtick);
    plot(Stick(:,1),Stick(:,2),'k')
end

%plot spokes representing constant declination
D0=[0:15:345];
for i=1:numel(D0)
    I1=[0 90]';
    D1=ones(size(I1)).*D0(i);
    XYZtick=ID2XYZ(deg2rad(I1),deg2rad(D1));
    Stick=pmag_splot(XYZtick);
    plot(Stick(:,1),Stick(:,2),'k')
end

set(gca,'box','off','visible','off')

%plot some random data
XYZ=randn(20,3);
XYZ=bsxfun(@rdivide,XYZ,sqrt(sum(XYZ.^2,2)));
Sxy=pmag_splot(XYZ);
plot(Sxy(:,1),Sxy(:,2),'or')
