function [U,MAD,score]=pmag_pca(x,y,z)

X=[x(:),y(:),z(:)]; %combine coordinates into single matrix

[V0,D0]=eigs(cov(X)); %leading eigenvectors & eigenvalues

%take leading eigenvector and set its orientation
U=sign(dot(mean(X)./norm(mean(X)),V0(:,1)))*V0(:,1); 

D0=diag(D0); %select eigenvalues
MAD=atan(sqrt((D0(2)+D0(3))./D0(1)))./pi.*180; %maximum angluar deviation (degrees)

score=bsxfun(@minus,X,mean(X))*U; %scores on leading PC


