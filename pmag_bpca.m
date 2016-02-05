function [MAP_I,MAP_D,theta95,MD0,Ip,Dp,VN,VE,VH, Emu] = pmag_bpca(x,y,z,p)

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
npts=21; % number of points along line
EZi=linspace(min(EZ)*1.1,max(EZ)*1.1,npts); %scores along line

 for i=1:npts
    [VNH,VNV,VEH,VEV,VHH,VHV]=pmag_dplane(EW,CW,Emu,Cmu,EZi(i),p);
    
    VN.H{1}(i,:)=VNH(1,:).*sqrt(Dn(1));
    VN.H{2}(i,:)=VNH(2,:).*sqrt(Dn(1));
    VN.V{1}(i,:)=VNV(1,:).*sqrt(Dn(1));
    VN.V{2}(i,:)=VNV(2,:).*sqrt(Dn(1));
 
    VE.H{1}(i,:)=VEH(1,:).*sqrt(Dn(1));
    VE.H{2}(i,:)=VEH(2,:).*sqrt(Dn(1));
    VE.V{1}(i,:)=VEV(1,:).*sqrt(Dn(1));
    VE.V{2}(i,:)=VEV(2,:).*sqrt(Dn(1));
    
    VH.H{1}(i,:)=VHH(1,:).*sqrt(Dn(1));
    VH.H{2}(i,:)=VHH(2,:).*sqrt(Dn(1));
    VH.V{1}(i,:)=VHV(1,:).*sqrt(Dn(1));
    VH.V{2}(i,:)=VHV(2,:).*sqrt(Dn(1));    
 end
 Emu = Emu.*sqrt(Dn(1));
