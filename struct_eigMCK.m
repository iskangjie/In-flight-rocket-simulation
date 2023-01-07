function [modes,freq,damping]=struct_eigMCK(M,C,K)


[ndof,tmp]=size(M);

A=[zeros(ndof,ndof) M;M C];

B=[-M zeros(ndof,ndof);zeros(ndof,ndof) K];

[eigvector,eigvalue]=eig(B,A);

freq=diag(imag(eigvalue)/(2*pi),0);
[freq,index]=sort(freq);
damp_factor=diag(real(eigvalue),0);
damp_factor=damp_factor(index);
damping=damp_factor./sqrt((2*pi*freq).^2+damp_factor.^2);

modes=eigvector(:,index);


