
%%GSTRPCA: Irregular tensor singular value decomposition for single-cell multi-omics
% data clustering
%%Input: Load irregular data types and construct irregular tensors
%        Parameter miu
%        Weight a
%%Output evaluation indicators:  ACC(Accuracy)
%                                ARI(Adjusted Rand Index)
%                                AMI(Adjusted Mutual Information)
%                                NMI(Normalized Mutual Information)
% We developed a weighted threshold model for the decomposition of irregular tensor data 
% by combining low-rank and sparsity constraints


clc;
clear;
miu=1e-4;


% init. variables  
 %load('10X_inhouse')%Loading data
 load('Specter')%Loading data
Tensor_X(:,:,1)=D1;
[m1,~]=size(D1);
[m2,~]=size(D2);
Tensor_X(m1+1:m2,:,1)=0;
Tensor_X(:,:,2)=D2;
Tensor_X(:,:,1)=sparse(Tensor_X(:,:,1));
Tensor_X(:,:,2)=sparse(Tensor_X(:,:,2));
yy1=zeros(20,20);

chu1=zeros(40,400);

kk=1;
j=1;
%iter
for p=0:0.05:0.95
    k=1;
      for a=0:0.05:0.95
          %Alternating iteration updates the main variables
[qq1,Tensor_epsilon1,Tensor_L1] = TSN1(miu,D1,D2,Tensor_X,p);
[m1,~]=size(qq1);
chu1(1:m1,kk)=qq1;
Tensor_L1=real(Tensor_L1);
Tensor_L1= max(Tensor_L1,0);
Tensor_L1(:,:,1)=sparse(Tensor_L1(:,:,1));
Tensor_L1(:,:,2)=sparse(Tensor_L1(:,:,2));


Y11=Tensor_L1(:,:,1)+a*Tensor_epsilon1(:,:,1);
Y12=Tensor_L1(:,:,2)+a*Tensor_epsilon1(:,:,2);


kk=kk+1;
%
Q1=[Y11(1:m1,:,1) 
    Y12];

%
for i=1:m1+m2
    SS1(i,:)=sum(Q1(i,:));
end
[row1,col1]=find(SS1<=1e-4);

Q1(row1,:)=[];

Q1(Q1<0)=0;

Q1=real(Q1);

%Calculate evaluation indicators

[ari,ami,nmi,accuracy] = SS(Q1);
yy1(k,j)=accuracy;
uu1(k,j)=ari;
ii1(k,j)=ami;
oo1(k,j)=nmi;


k=k+1;
      end
      j=j+1;
end
