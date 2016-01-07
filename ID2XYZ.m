function xyz=ID2XYZ(I,D)

X=cos(D).*cos(I);
Y=sin(D).*cos(I);
Z=sin(I);
xyz=[X Y Z];