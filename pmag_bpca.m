function [MAP_I,MAP_D,theta95,MD0,Ip,Dp,VN,VE,VH,Emu] = pmag_bpca(x,y,z,p)

q=1;
X=[x(:)';y(:)';z(:)'];
[m,n] = size(X);

%% normalize data to approximatly correct scale
Dn=eigs(cov(X'));
X=X./sqrt(Dn(1));

%% BPCA
[V0,D0]=eigs(cov(X'));
V0(:,1)=sign(dot(mean(X,2)./norm(mean(X,2)),V0(:,1)))*V0(:,1);
D0=diag(D0);

a0 = 1e-3;
b0 = 1e-3;
c0 = 1e-3;
d0 = 1e-3;
lambda=1e-3;
maxIter = 1e4;

a = a0+m/2;
c = c0+m*n/2;
Ealpha = 1e-1;

Ebeta = 1./mean(D0(2:3));
EW = V0(:,1:q)'.*sqrt(D0(1)-mean(D0(2:3)));
EWo = bsxfun(@minus,EW,mean(EW,2));
EWW = EWo*EWo'/m+EW*EW';
I = eye(q);

Emu = mean(X,2);
Xo = bsxfun(@minus, X, Emu);
s = dot(Xo(:),Xo(:));

for iter = 1:maxIter
    % calculate z
    CZ = inv(I+Ebeta*EWW);
    EZ = Ebeta*CZ*EW*Xo;
    EZZ = n*CZ+EZ*EZ';
    
    % calculate w
    A = diag(Ealpha);
    CW = inv(A+Ebeta*EZZ);
    EW = Ebeta*CW*EZ*Xo';
    EWW = m*CW+EW*EW';
    
    % calculate alpha
    b = b0+diag(EWW)/2;
    Ealpha = a./b;
    
    % calculate beta
    WZ = EW'*EZ;
    d = d0+(s-2*dot(Xo(:),WZ(:))+dot(EWW(:),EZZ(:)))/2;
    Ebeta = c/d;
    
    % calculate mu
    Cmu = 1./(lambda+n*Ebeta);
    Emu = Ebeta*Cmu*sum(X-WZ,2);
    Xo = bsxfun(@minus, X, Emu);
    s = dot(Xo(:),Xo(:));
end

%% calculate MAP direction and statistics
[MAP_I,MAP_D]=XYZ2ID(EW);
MAP_I=rad2deg(MAP_I);
MAP_D=rad2deg(MAP_D);

%% calculate theta95;
theta95=rad2deg(agauss_inv(EW,CW,p));

%% test the origin

%%%Emu is expected mean(rescale at last)
for i=1:1e4
    W1=mvnrnd(EW,eye(3)*CW);
    M1=mvnrnd(Emu,eye(3)*Cmu);
    I1(i,:)=-M1*pinv(W1)*W1+M1;
end
Zbar=mean(I1);
ZS=cov(Zbar);
MD0=-Zbar*pinv(ZS)*-Zbar';

%% marginalize directions
t=mvnrnd(EW,eye(3)*CW,1e4);
[It,Dt]=XYZ2ID(t);
Ip=rad2deg(prctile(It,[(1-p)./2 1-(1-p)./2]*100));
Dp=rad2deg(prctile(Dt,[(1-p)./2 1-(1-p)./2]*100));

%% %% calculate XY HDI
 npts=51; % number of points along line
% 
xerr12=NaN(npts,1);
xerr21=NaN(npts,1);
yerr12=NaN(npts,1);
yerr21=NaN(npts,1);
xerr13=NaN(npts,1);
xerr31=NaN(npts,1);
zerr13=NaN(npts,1);
zerr31=NaN(npts,1);
% 
 EZi=linspace(min(EZ)*1.1,max(EZ)*1.1,npts); %scores along line
for i=1:npts
    bar=bsxfun(@plus,EZi([1 i end])'*EW,Emu');
    [~,~,Dlim]=pmag_dplane(EW,CW,Emu,Cmu,EZi(i),p);
    
    pbar(1)=(bar(end,2)-bar(1,2))./(bar(end,1)-bar(1,1));
    pbar(2)=bar(1,2)-pbar(1)*bar(1,1);
        
    x0=bar(2,1);
    y0=bar(2,2);
    
    pnorm=[-1./pbar(1) y0+x0./pbar(1)];
    
    xerr12(i,1)=x0+sqrt(Dlim.^2./(1+pnorm(1).^2));
    yerr12(i,1)=pnorm(1)*xerr12(i,1)+pnorm(2);
    xerr21(i,1)=x0-sqrt(Dlim.^2./(1+pnorm(1).^2));
    yerr21(i,1)=pnorm(1)*xerr21(i,1)+pnorm(2);
    
    % calculate XZ HDI
    pbar(1)=(bar(end,3)-bar(1,3))./(bar(end,1)-bar(1,1));
    pbar(2)=bar(1,3)-pbar(1)*bar(1,1);
    
    x0=bar(2,1);
    z0=bar(2,3);
    pnorm=[-1./pbar(1) z0+x0./pbar(1)];
    xerr13(i,1)=x0+sqrt(Dlim.^2./(1+pnorm(1).^2));
    zerr13(i,1)=pnorm(1).*xerr13(i,1)+pnorm(2);
    xerr31(i,1)=x0-sqrt(Dlim.^2./(1+pnorm(1).^2));
    zerr31(i,1)=pnorm(1).*xerr31(i,1)+pnorm(2);
 end

x12=[xerr12;flipud(xerr21);xerr12(1)];
y12=[yerr12;flipud(yerr21);yerr12(1)];
x13=[xerr13;flipud(xerr31);xerr13(1)];
z13=[zerr13;flipud(zerr31);zerr13(1)];

%% rescale to original data magnitude
bar=bar.*sqrt(Dn(1));
x12=x12.*sqrt(Dn(1));
y12=y12.*sqrt(Dn(1));
x13=x13.*sqrt(Dn(1));
z13=z13.*sqrt(Dn(1));

xerr12=xerr12.*sqrt(Dn(1));
xerr21=xerr21.*sqrt(Dn(1));
yerr12=yerr12.*sqrt(Dn(1));
yerr21=yerr21.*sqrt(Dn(1));

xerr13=xerr13.*sqrt(Dn(1));
xerr31=xerr31.*sqrt(Dn(1));
zerr13=zerr13.*sqrt(Dn(1));
zerr31=zerr31.*sqrt(Dn(1));

%% package into VN outputs
VN.H{1}=[yerr12,xerr12];
VN.H{2}=[yerr21,xerr21];
VN.V{1}=[xerr13,-zerr13];
VN.V{2}=[xerr31,-zerr31];

%% package into VE outputs
VE=VN;

%% package in  VH outputs
VH=VN;


Emu = Emu.*sqrt(Dn(1));
