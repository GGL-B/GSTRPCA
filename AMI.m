function [AMI_]=AMI(true_mem,mem)
%Program for calculating the Adjusted Mutual Information (AMI) between
%two clusterings.用于计算两个聚类之间的调整后互信息（AMI）的程序。
%
% This code is a modified version of Nguyen Xuan Vinh's one available online.
%
% The modification includes the computation of AMI for contingency tables
% with empty rows or columns.修改包括计算具有空行或空列的列联表的AMI。
%
%--------------------------------------------------------------------------
%**Input: a contingency table T，输入：应急表T或者两个向量中两个聚类的聚类标记
%   OR
%        cluster label of the two clusterings in two vectors
%        eg: true_mem=[1 2 4 1 3 5]
%                 mem=[2 1 3 1 4 5]
%        Cluster labels are coded using positive integer. 集群标签使用正整数编码。
%**Output: AMI: adjusted mutual information  (AMI normalized by Sqrt(HA,HB))，输出：AMI：调整后的相互信息（通过Sqrt（HA，HB）归一化的AMI）
%
%**Note: In a prevous published version, if you observed strange AMI results, eg. AMI>>1, 
%then it's likely that in these cases the expected MI was incorrectly calculated (the EMI is the sum
%of many tiny elements, each falling out the precision range of the computer).
%However, you'll likely see that in those cases, the upper bound for the EMI will be very
%tiny, and hence the AMI -> NMI (see [3]). It is recommended setting AMI=NMI in
%these cases, which is implemented in this version.
%在以前发布的版本中，如果你观察到奇怪的AMI结果，例如AMI>>1，那么在这些情况下，预期的MI很可能计算错误（EMI是许多微小元素的总和，每个都超出了计算机的精度范围）。
%然而，您可能会看到，在这些情况下，EMI的上限将非常小，因此AMI->NMI（参见[3]）。在这些情况下，建议设置AMI=NMI，这在本版本中实现。
%References: 
% [1] 'A Novel Approach for Automatic Number of Clusters Detection based on Consensus Clustering', 
%       N.X. Vinh, and Epps, J., in Procs. IEEE Int. Conf. on 
%       Bioinformatics and Bioengineering (Taipei, Taiwan), 2009.
% [2] 'Information Theoretic Measures for Clusterings Comparison: Is a
%	    Correction for Chance Necessary?', N.X. Vinh, Epps, J. and Bailey, J.,
%	    in Procs. the 26th International Conference on Machine Learning (ICML'09)
% [3] 'Information Theoretic Measures for Clusterings Comparison: Variants, Properties, 
%       Normalization and Correction for Chance', N.X. Vinh, Epps, J. and
%       Bailey, J., Journal of Machine Learning Research, 11(Oct), pages
%       2837-2854, 2010


if nargin==1
    T=true_mem; %contingency table pre-supplied，预先提供
elseif nargin==2
    %build the contingency table from membership arrays，从成员数组构建列联表
    r=max(true_mem);
    c=max(mem);

    %identify & removing the missing labels，识别并删除丢失的标签
    list_t=ismember(1:r,true_mem);%ismember这个函数主要是看矩阵A中的数据是不是矩阵B中的成员，是的话返回一个包含逻辑1(ture)的数组，不是返回0；
    list_m=ismember(1:c,mem);
    T=Contingency(true_mem,mem);
    T=T(list_t,list_m);
end

[r c]=size(T);
if (c == 1 || r == 1)
 error('Clusterings should have at least 2 clusters')
 return
end

N = sum(sum(T)); % total number of records，记录总数

% update the true dimensions，更新真实尺寸
a=sum(T,2)';
b=sum(T);

%calculating the Entropies，计算熵
Ha=-(a(a ~= 0)/N)*log2(a(a ~= 0)/N)'; 
Hb=-(b(b ~= 0)/N)*log2(b(b ~= 0)/N)';

%calculate the MI (unadjusted)，计算MI（未调整）
MI=0;
for i=1:r
    for j=1:c
        if T(i,j)>0 MI=MI+T(i,j)*log2(T(i,j)*N/(a(i)*b(j)));end;
    end
end
MI=MI/N;

%-------------correcting for agreement by chance，偶然纠正---------------------------
AB=a'*b;
bound=zeros(r,c);

E3=(AB/N^2).*log2(AB/N^2);
E3(isnan(E3)) = 0; % substitute 0log0=NaN with 0s

EPLNP=zeros(r,c);
log2Nij=log2([1:min(max(a),max(b))]/N);
for i=1:r
    for j=1:c
        nij=max(1,a(i)+b(j)-N);
        X=sort([nij N-a(i)-b(j)+nij]);
        if N-b(j)>X(2)
            nom=[[a(i)-nij+1:a(i)] [b(j)-nij+1:b(j)] [X(2)+1:N-b(j)]];
            dem=[[N-a(i)+1:N] [1:X(1)]];
        else
            nom=[[a(i)-nij+1:a(i)] [b(j)-nij+1:b(j)]];       
            dem=[[N-a(i)+1:N] [N-b(j)+1:X(2)] [1:X(1)]];
        end
        p1=prod(nom./dem)/N;                
        for nij=max(1,a(i)+b(j)-N):1:min(a(i), b(j))
            EPLNP(i,j)=EPLNP(i,j)+nij*log2Nij(nij)*p1;            
            p1=p1*(a(i)-nij)*(b(j)-nij)/(nij+1)/(N-a(i)-b(j)+nij+1);  
        end
        CC=N*(a(i)-1)*(b(j)-1)/a(i)/b(j)/(N-1)+N/a(i)/b(j);
        bound(i,j)=a(i)*b(j)/N^2*log2(CC);         
    end
end

EMI_bound=sum(sum(bound));
EMI_bound_2=log2(r*c/N+(N-r)*(N-c)/(N*(N-1)));
EMI=sum(sum(EPLNP-E3));

AMI_=(MI-EMI)/(sqrt(Ha*Hb)-EMI);
NMI=MI/sqrt(Ha*Hb);


%If expected mutual information negligible, use NMI.如果预期互信息可忽略不计，则使用NMI。
if abs(EMI)>EMI_bound
    fprintf('The EMI is small: EMI < %f, setting AMI=NMI',EMI_bound);
    AMI_=NMI;
end;

%---------------------auxiliary functions，辅助功能---------------------
function Cont=Contingency(Mem1,Mem2)

if nargin < 2 || min(size(Mem1)) > 1 || min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end