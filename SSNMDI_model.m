function [U_final,Z_final,B_final,F_final, S_final] = SSNMDI_model(X, A, lambda, k1, k2, W, options, num)

differror = 1e-5; %
if isfield(options,'error') 
    differror = options.error;
end

maxIter = []; 
if isfield(options, 'maxIter')
    maxIter = options.maxIter;
end

nRepeat =1;
if isfield(options,'nRepeat')
    nRepeat = options.nRepeat;
end

minIterOrig = 30; 
if isfield(options,'minIter')
    minIterOrig = options.minIter;
end
minIter = minIterOrig-1;

meanFitRatio = 0.1;
if isfield(options,'meanFitRatio')
    meanFitRatio = options.meanFitRatio;
end

alpha = 1;
if isfield(options,'alpha')
    alpha = options.alpha;
end

if min(min(X)) < 0
    error('Input should be nonnegative!');
end

[m,n]=size(X);
U0= abs(rand(m,k1));
Z0 = abs(rand(n+k2-num,k1));%U0,Z0
[mFea,nSmp]=size((A * Z0)');%mFea=k1,nSmp=n
%%%%%%=========== make S0===========
window=(X==0); 
dim=size(X);
rand('seed',23);
S0=rand(dim).*window;
%%%%%%=========== Weight matrix setting===========

if isfield(options,'weight') && strcmpi(options.weight,'NCW') %strcmpi
    feaSum = full(sum(X,2)); %X:n*m,feaSum:n*1)
    D_half = (X'*feaSum).^.5;%m*n*n*1=m*1
    for i = 1:nSmp
        X(:,i) = X(:,i)/D_half(i);
    end
end

if isfield(options,'alpha_nSmp') && options.alpha_nSmp
    alpha = alpha*nSmp;    
end

W = alpha*W;% Weight matrix

DCol = full(sum(W,2));% Sum of row elements constitutes column vector DCol  
D = spdiags(DCol,0,speye(size(W,1)));% Compose Diagonal D 
L = D - W;
if isfield(options,'NormW') && options.NormW 
    D_mhalf = DCol.^-.5;

    tmpD_mhalf = repmat(D_mhalf,1,nSmp);
    L = (tmpD_mhalf.*L).*tmpD_mhalf';
    clear D_mhalf tmpD_mhalf;

    L = max(L, L');
end


%%%%%%%===========initialization================
    
selectInit = 1;
if ~exist('U','var')
    B0 = abs(rand(k1,k2));
    F0 = abs(rand(k2,n));
else
    nRepeat = 1;
end
%%%%%%===========Parameter settings===========
Norm = 2;
NormV = 1;
[U0,Z0] = NormalizeUZ(U0, Z0, NormV, Norm);
[B0,F0] = NormalizeUV(B0, F0', NormV, Norm);F0=F0';

Uk=U0;Zk=Z0;
Bk=B0;Fk=F0;
Sk = S0;
Ek=Zk;
[mFea1,nSmp1]=size(Zk);
Tk= zeros(mFea1,k1);

iter = 0; 
converged = 0;    
maxIter=100;  
tol1=1e-5;tol2=1e-5;
tryNo=0;
w1=1;w2=1;

%%%%%%%===========Update variables U,Z,B,F,S by iteration================

  while ~converged  && iter < maxIter   
    iter = iter + 1;
        derta =5e+1;

        alpha=norm(X,1)/norm((A*Zk)','fro');beta=norm((A*Zk)',1)/norm(Fk,'fro');%%%%%%%Regularization parameter 
        E=eye(k1);
        %%%========Update variables U,Z,B==========
        Ukl=Uk .* ( ((X+Sk)*A*Zk) ./ (Uk*Zk'*A'*A*Zk+alpha*Uk) );
        Zkl=Zk .* ( (A'*(X+Sk)'*Uk+A'*Fk'*Bk') ./ (A'*A*Zk*(Uk'*Uk+E)) );
%         I=eye(k1);
%         VV1=A'*(X+Sk)'*Uk+A'*Fk'*Bk'+derta*Ek+Tk;
%         VV2=A'*A*Zk*(Uk'*Uk+E)+Zk*I+derta*Zk;
%         Zkl=Zk.*(VV1./VV2);
        Ekl=soft(Zkl-Tk/derta,alpha/derta);
        Tkl=Tk+1.618*derta*(Ekl-Zkl);
        Bkl=Bk .* ( (Zkl'*A'*Fk') ./ (Bk*Fk*Fk') );
        for i=1:size(Bkl,1) 
            for j=1:size(Bkl,2)
                if Bkl(i,j)<0
                   Bkl(i,j)=0;
                else
                    Bkl(i,j)=Bkl(i,j);
                end
            end
        end
        %%%========Update variables F==========
        FF1=Bkl'*Zkl'*A'+beta*Fk*W;
        FF2=Bkl'*Bkl*Fk+beta*Fk*D;
        Fkl=Fk.*(FF1)./(FF2);
        for i=1:size(Fkl,1)
            for j=1:size(Fkl,2)
                if Fkl(i,j)<0
                   Fkl(i,j)=0;
                else
                    Fkl(i,j)=Fkl(i,j);
                end
            end
        end
        %%%========Update variables S==========
        temp=-X+Ukl*Zkl'*A';
        miu=diag(sum(X)/median(sum(X)));
        save temp temp;
        save miu miu;
        Skl=SoftThreshold(temp,ones(dim)*(lambda*miu));
        Skl=Skl.*window;
%         [DTkl,Vkl] = NormalizeUV(DTkl, Vkl', NormV, Norm);Vkl=Vkl';
        [Bkl,Fkl] = NormalizeUV(Bkl, Fkl', NormV, Norm);Fkl=Fkl';
    
    Uwk = Uk;
    Zwk = Zk;
    Bwk = Bk;
    Fwk = Fk;
    Swk = Sk; 
    
    Uk = Ukl;
    Zk = Zkl;
    Bk = Bkl;
    Fk = Fkl;
    Sk = Skl;
%%%%%%%%%% Error 
  Er1(iter,:)=abs(mean((X - Uk*Zk'*A')))./norm((Uk*Zk'*A'),'fro');
  Er2(iter,:)=abs(mean((Zk'*A' - Bk*Fk)))./norm( Bk*Fk,'fro');
  er1(iter)=mean(Er1(iter,:));
  er2(iter)=mean(Er2(iter,:));
  er=er1+er2;
  

    temp = max ([norm(Ukl-Uwk,'fro'),norm(Zkl-Zwk,'fro'),norm(Bkl-Bwk,'fro'),norm(Fkl-Fwk,'fro')]);
%         temp = muu*temp/norm(V,2);
%     temp =temp/max([norm(X,'fro')]);
temp =temp/norm(X,'fro');
%         temp = max([(sqrt(L)*norm(ZK-Zkm1,'fro')),norm(WK-Wkm1,'fro'),norm(EK-Ekm1,'fro')]);
%         temp = muu*temp/norm(Y,'fro');
        
    %%%%%%%%%%%%%%%%%%
    temp1 = max(norm( (X - Uk*Zk'*A'),'fro'),norm( (Zk'*A' - Bk*Fk),'fro'))/max(norm( Bk*Fk,'fro'),norm( Uk*Zk'*A','fro'));
    if temp1 < tol1 && temp < tol2
    converged = 1;
    end
     disp(['temp1 ',num2str(temp1)]);
     disp([' µü´ú´ÎÊý ' num2str(iter) ' temp1 ' num2str(temp1) ' temp ' num2str(temp)]);
    t1(iter)=temp1;
    t2(iter)=temp;
%     elapse = cputime - tmp_T;

  end
  
        U_final = Ukl; 
        Z_final = Zkl;
        B_final = Bkl;
        F_final = Fkl; 
        S_final = Skl;
        
        [U_final,Z_final] = NormalizeUZ(U_final, Z_final, NormV, Norm);
        [B_final,F_final] = NormalizeUV(B_final, F_final', NormV, Norm);F_final=F_final';
%         t=1:iter;
%         figure
%         plot(t,er,'r-'),xlabel('Iteration times');ylabel('Relative Error');
%         hold on;
end

 function[y] = soft( x, T )
   if sum( abs(T(:)) )==0
        y = x;
   else
        y = max( abs(x) - T, 0);
        y = sign(x).*y;
    end
 end    