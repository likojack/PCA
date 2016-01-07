function XY=pmag_splot(XYZ)

[I,D]=XYZ2ID(XYZ);
D=D+pi./2;
XYZ=ID2XYZ(I,D);
XY=XYZ(:,1:2);
XY=bsxfun(@rdivide,XY,sqrt(sum(XY.^2,2)));
XY(:,1)=-XY(:,1);

L=(pi./2-asin(abs(XYZ(:,3))))./(pi./2);
XY=bsxfun(@times,XY,L);