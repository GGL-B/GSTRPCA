function [Tensor_L] =  Imsvd1(Y,rho)


epsilon=1e-8;
p=0;
[n1,n2,n3] = size(Y);
Tensor_L = zeros(n1,n2,n3);

%%Fourier transform
Y = fft(Y,[],3);

        
% first frontal slice
%SVD
[U,S,V] = svd(Y(:,:,1),'econ');
S = diag(S);
r = length(find(S>rho));
if r>=1
    miu3=diag(S(1:r))-1./((S(1:r)-epsilon).^(1-p));
 C=diag(max(miu3,0));
    Tensor_L(:,:,1) = U(:,1:r)*diag(C)*V(:,1:r)';

end

halfn3 = round(n3/2);

% second frontal slice
%SVD
if mod(n3,2) == 0
    i = halfn3+1;
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    r = length(find(S>rho));
    if r>=1
                   miu3=diag(S(1:r))-1./((S(1:r)-epsilon).^(1-p));
       C=diag(max(miu3,0));
        Tensor_L(:,:,i) = U(:,1:r)*diag(C)*V(:,1:r)';
  
    end
end
Tensor_L= ifft(Tensor_L,[],3);
