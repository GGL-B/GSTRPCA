function [S,B1,B2,V,U1,V1,S1,eps,eps2,sigma,c,U,A,Z]=TwoSVD(D1,D2,KK)

tsize = size(D1);
Z    = tsize(2);

%Construct S
A(:,:,1)=D1.'*D1;
A(:,:,2)=D2.'*D2;

    S=1/2*(A(:,:,1)*pinv(A(:,:,2))+A(:,:,2)*pinv(A(:,:,1)));
    
    [V,~]=eig(S);
    B1=V^(-1)*D1.';
    B2=V^(-1)*D2.';
    B1=B1.';
    B2=B2.';


    %caculate sigma£¬caculate U
   for i=1:Z
       if KK==1
    c(1,i)=norm(B1(:,i));
       sigma=diag(c(1,:)); 
       else
        c(1,i)=norm(B2(:,i));    
           sigma=diag(c(1,:)); 
       end
   end
  
   if KK==1
   U=B1*pinv(sigma);
     eps=U*sigma*V.'-D1;
      [U1,S1,V1]=svd(D1);
   else
          U=B2*pinv(sigma); 
             eps=U*sigma*V.'-D2;
               [U1,S1,V1]=svd(D2);
   end
  
 

eps2=U1*S1*V1.'-U*sigma*V.';

end