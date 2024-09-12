function[A] = labelA(label, c, k2)
n = size(label, 1);%该语句返回label的行数
% N=randi(n,1,l);
num = length(c);
C = zeros(num, k2);
for i=1:num
    C(i, c(i))=1;
end
I = eye(n-num);
A = zeros(n, n + k2 - num);
for i=1:num
    for j=1:k2
        A(i, j) = C(i, j);
    end
end
for i2=1:n-num
        A(num + i2, k2 + i2) = I(i2, i2);
end