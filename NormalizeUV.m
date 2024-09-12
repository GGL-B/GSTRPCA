function [U, V] = NormalizeUV(U, V, NormV, Norm)
%NormV，V的2范数
    nSmp = size(V,1);%该语句返回V的行数，2是列数
    mFea = size(U,1);%该语句返回U的行数
    if Norm == 2
        if NormV
            norms = sqrt(sum(V.^2,1));%sum(A,1or2),1表示每一列进行求和，2表示每一行进行求和；
            norms = max(norms,1e-10);
            V = V./repmat(norms,nSmp,1);
            U = U.*repmat(norms,mFea,1);
            %norms(a,1,1) a的列求l1范数，结果是行向量
            %norms(a,1,2) a的行求l1范数，结果是列向量
            %norms(a,1,3) 每个数求l1范数
            %norms(a,1) 没有s，就是矩阵的范数，最大列和
        else
            norms = sqrt(sum(U.^2,1));
            norms = max(norms,1e-10);
            U = U./repmat(norms,mFea,1);
            V = V.*repmat(norms,nSmp,1);
        end
    else
        if NormV
            norms = sum(abs(V),1);
            norms = max(norms,1e-10);%原矩阵负元素变为1e-10
            V = V./repmat(norms,nSmp,1);
            U = U.*repmat(norms,mFea,1);%repmat：复制和平铺矩阵，将norms复制平铺为以norms为单位的mFea*1阶矩阵
        else
            norms = sum(abs(U),1);
            norms = max(norms,1e-10);%原矩阵负元素变为1e-10
            U = U./repmat(norms,mFea,1);
            V = V.*repmat(norms,nSmp,1);
        end
    end

end    