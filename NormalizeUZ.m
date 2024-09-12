function [U, Z] = NormalizeUV(U, Z, NormV, Norm)
    nSmp = size(Z,1);
    mFea = size(U,1);
    if Norm == 2
        if NormV
            norms = sqrt(sum(Z.^2,1));
            norms = max(norms,1e-10);%原矩阵负元素变为1e-10
            Z = Z./repmat(norms,nSmp,1);%repmat是复制和平铺矩阵
            U = U.*repmat(norms,mFea,1);
        else
            norms = sqrt(sum(U.^2,1));
            norms = max(norms,1e-10);
            U = U./repmat(norms,mFea,1);
            Z = Z.*repmat(norms,nSmp,1);
        end
    else
        if NormV
            norms = sum(abs(Z),1);
            norms = max(norms,1e-10);
            Z = Z./repmat(norms,nSmp,1);
            U = U.*repmat(norms,mFea,1);
        else
            norms = sum(abs(U),1);
            norms = max(norms,1e-10);
            U = U./repmat(norms,mFea,1);
            Z = Z.*repmat(norms,nSmp,1);
        end
    end

end    