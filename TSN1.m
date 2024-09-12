function [qq1,Tensor_epsilon,Tensor_L,Tensor_Y] = TSN1(miu,D1,D2,Tensor_X,p)
% TSN is the main running program of GSTRPCA,
%%Input: Load irregular data types and Construct matrices D1 and D2 
%        Parameter miu
%        Weight p
%        Irregular tensor X
%%Output: Irregular low-rank tensor L
%         Irregular sparse tensor
%         Lagrange Multiplier Tensor Y


[m1,~]=size(D1);
[m2,~]=size(D2);
tsize = size(Tensor_X);
N    = tsize(1);
M    = tsize(2);
K    = tsize(3);
lambda=10/sqrt(max(N,M)*K);
Tensor_Y=zeros(tsize);
Tensor_epsilon=zeros(tsize);
Tensor_L=zeros(tsize);
pro=1.1;
tol=1e-3;
miumax=1e+5;
maxiter=500;
display  = 1;
k=1;
for iter=1:maxiter
     fprintf('\n*****************************iter: %d ******************************\n', iter');
       L_pre=Tensor_L;
       epsilon_pre=Tensor_epsilon;
      bbb=real(Tensor_X-Tensor_epsilon-Tensor_Y/miu);     
        
  %Irregular low-rank tensor L iteration

[U,V,sigma,Tensor_L] =t_SVT1(bbb,1/miu,p,m1,m2);

%Irregular sparse tensor epsilon iteration
  Tensor_epsilon = prox_l1(Tensor_X-Tensor_L-Tensor_Y/miu,lambda/miu);

    dY = Tensor_L+Tensor_epsilon-Tensor_X;
    chgL = max(abs(L_pre(:)-Tensor_L(:)));
    chgS = max(abs(epsilon_pre(:)-Tensor_epsilon(:)));
    chg = max([ chgL chgS max(abs(dY(:))) ]);
    err = norm(dY(:));
    
%Iterative convergence
  if  display
        fprintf('miu:%4.4e,chg:%4.4e,chgL:%4.4e, chgS:%4.4e,err: %4.4e\n', miu,chg,chgL,chgS,err);
  end
   if err < tol 
          disp(' !!!stopped by termination rule!!! ');  break;
   end
   
%Lagrange Multiplier Tensor Y iteration
    Tensor_Y=Tensor_Y+miu*(dY);
    Tensor_Y=real(Tensor_Y);
    
% Step size miu iteration
    miu=min(pro*miu,miumax);  
    
    if mod(iter,1)==0
        qq1(k,:)=err;
        k=k+1;
    end
    
end


end

