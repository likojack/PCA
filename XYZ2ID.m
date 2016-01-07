function [I,D]=XYZ2ID(xyz)

xyz=bsxfun(@rdivide,xyz,sqrt(sum(xyz.^2,2)));
I=asin(xyz(:,end));
D=atan2(xyz(:,2),xyz(:,1));