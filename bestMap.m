function [newL2] = bestMap(L1,L2)
%bestmap: permute labels of L2 to match L1 as good as possible  bestmap：排列L2的标签以尽可能好地匹配L1
%   [newL2] = bestMap(L1,L2);
%
%   version 2.0 --May/2007
%   version 1.0 --November/2003
%
%   Written by Deng Cai (dengcai AT gmail.com)


%===========    

L1 = L1(:);%将L1列表按列展开，是一个列向量
L2 = L2(:);;%将L2列表按列展开，是一个列向量
if size(L1) ~= size(L2)
    error('size(L1) must == size(L2)');
end

Label1 = unique(L1);%返回与 L1 中相同的数据，但是不包含重复项
nClass1 = length(Label1);
Label2 = unique(L2);%返回与 L2 中相同的数据，但是不包含重复项
nClass2 = length(Label2);

nClass = max(nClass1,nClass2);
G = zeros(nClass);%返回一个nClass x nClass的零矩阵
for i=1:nClass1
	for j=1:nClass2
		G(i,j) = length(find(L1 == Label1(i) & L2 == Label2(j)));
	end
end

[c,t] = hungarian(-G);
newL2 = zeros(size(L2));
for i=1:nClass2
    newL2(L2 == Label2(i)) = Label1(c(i));
end

