function [modes,freq]=struct_eigMK(M,K)
% [eigvector,eigvalue]=eig(K,M,'qz');
[eigvector,eigvalue]=eig(K,M);
freq=sqrt(diag(eigvalue,0))/(2*pi);
modes=eigvector;

