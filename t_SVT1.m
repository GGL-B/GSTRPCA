function [U,V,sigma,Tensor_L] =t_SVT1(bbb,rho,p,m1,m2)
%%no-fill
C1=bbb(1:m1,:,1);
 C2=bbb(:,:,2);
epsilon=1e-8;
 tsize = size(bbb);
K    = tsize(3);
Tensor_L=zeros(tsize);

  %%%first frontal slice
halfn3 = round(K/2);
for i = 1 : halfn3
    
[~,~,~,V,~,~,~,~,~,sigma,~,U,~,~]=TwoSVD(C1,C2,i);
  SS = diag(sigma);
  r = length(find(SS>rho));
if r>=1
miu3=diag(SS(1:r))-1./((SS(1:r)-epsilon).^(1-p));
 C=diag(max(miu3,0));
  U(m1+1:m2,:,1)=0;
  Tensor_L(:,:,i)=U(:,1:r)*diag(C)*V(:,1:r).';
end
end

% second frontal slice
if mod(K,2) == 0
    i = halfn3+1;

[~,~,~,V,~,~,~,~,~,sigma,~,U,~,~]=TwoSVD(C1,C2,i);

    SS = diag(sigma);
    r = length(find(SS>rho));
    if r>=1
       miu3=diag(SS(1:r))-1./((SS(1:r)-epsilon).^(1-p));
       C=diag(max(miu3,0));

 Tensor_L(:,:,i)=U(:,1:r)*diag(C)*V(:,1:r).';
    end
end

