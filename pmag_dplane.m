function [XYZ0,XYZi,D95]=pmag_dplane(EW,CW,Emu,Cmu,Z,p,iter)

if nargin==6
    iter=1e4;
end

%point on plane
XYZ0=Z*EW+Emu';

%normal vector to the plane is simply EW
n=EW; %[A,B,C]
D=dot(n,XYZ0);

mui=mvnrnd(Emu,eye(3)*Cmu,iter);
Wi=mvnrnd(EW,eye(3)*CW,iter);

for i=1:iter
    Nv=dot(n,Wi(i,:));
    Nr0=dot(n,mui(i,:));
    Zhat=(D-Nr0)./Nv;
    XYZi(i,:)=Zhat*Wi(i,:)+mui(i,:); %intersecting point
    Dmle(i,1)=norm(XYZi(i,:)-XYZ0);
end

D95=prctile(Dmle,p*100);
