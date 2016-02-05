function [VNH,VNV,VEH,VEV,VHH,VHV]=pmag_dplane(EW,CW,Emu,Cmu,Z,p,iter)

if nargin==6
    iter=5e3;
end

%point on plane
XYZ0=Z*EW+Emu';

%normal vector to the plane is simply EW
n=EW; %[A,B,C]
D=dot(n,XYZ0);

mui=mvnrnd(Emu,eye(3)*Cmu,iter);
Wi=mvnrnd(EW,eye(3)*CW,iter);

XYZi=zeros(iter,3);
for i=1:iter
    Nv=dot(n,Wi(i,:));
    Nr0=dot(n,mui(i,:));
    Zhat=(D-Nr0)./Nv;
    XYZi(i,:)=Zhat*Wi(i,:)+mui(i,:); %intersecting point
end


%% calculate limits on V vs N projection
% Horizontal (Y,X)
C=EW([2 1]).*[1 -1];
C=fliplr(C);
H=XYZi(:,[2 1]);
M=XYZ0([2 1]);
S=bsxfun(@minus,H,M)*pinv(C')';
pc=norminv(1-(1-p)./2);
VNH=bsxfun(@plus,[-pc*std(S);pc*std(S)]*C,M);

%vertical (X,-Z)
C=EW([1 3]);
C=fliplr(C);
V=XYZi(:,[1 3]);
V(:,2)=-V(:,2);
M=XYZ0([1 3]);
M(:,2)=-M(:,2);
S=bsxfun(@minus,V,M)*pinv(C')';
VNV=bsxfun(@plus,[-pc*std(S);pc*std(S)]*C,M);

%% calculate limits on V vs E projection
% Horizontal (Y,X)
VEH=VNH;

%vertical (Y,-Z)
C=EW([2 3]);
C=fliplr(C);
V=XYZi(:,[2 3]);
V(:,2)=-V(:,2);
M=XYZ0([2 3]);
M(:,2)=-M(:,2);
S=bsxfun(@minus,V,M)*pinv(C')';
VEV=bsxfun(@plus,[-pc*std(S);pc*std(S)]*C,M);

%% calculate limits on V vs H projection
% Horizontal (Y,X)
VHH=VNH;

%vertical (sqrt(X^2+Y^2),-Z)
C=[sqrt(EW(1).^2+EW(2).^2) EW(3)];
C=fliplr(C);
V(:,1)=sqrt(sum(XYZi(:,1:2).^2,2));
V(:,2)=-XYZi(:,3);
M(:,1)=sqrt(sum(XYZ0(:,1:2).^2,2));
M(:,2)=-XYZ0(:,3);
S=bsxfun(@minus,V,M)*pinv(C')';
VHV=bsxfun(@plus,[-pc*std(S);pc*std(S)]*C,M);
