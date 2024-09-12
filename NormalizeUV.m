function [U, V] = NormalizeUV(U, V, NormV, Norm)
%NormV��V��2����
    nSmp = size(V,1);%����䷵��V��������2������
    mFea = size(U,1);%����䷵��U������
    if Norm == 2
        if NormV
            norms = sqrt(sum(V.^2,1));%sum(A,1or2),1��ʾÿһ�н�����ͣ�2��ʾÿһ�н�����ͣ�
            norms = max(norms,1e-10);
            V = V./repmat(norms,nSmp,1);
            U = U.*repmat(norms,mFea,1);
            %norms(a,1,1) a������l1�����������������
            %norms(a,1,2) a������l1�����������������
            %norms(a,1,3) ÿ������l1����
            %norms(a,1) û��s�����Ǿ���ķ���������к�
        else
            norms = sqrt(sum(U.^2,1));
            norms = max(norms,1e-10);
            U = U./repmat(norms,mFea,1);
            V = V.*repmat(norms,nSmp,1);
        end
    else
        if NormV
            norms = sum(abs(V),1);
            norms = max(norms,1e-10);%ԭ����Ԫ�ر�Ϊ1e-10
            V = V./repmat(norms,nSmp,1);
            U = U.*repmat(norms,mFea,1);%repmat�����ƺ�ƽ�̾��󣬽�norms����ƽ��Ϊ��normsΪ��λ��mFea*1�׾���
        else
            norms = sum(abs(U),1);
            norms = max(norms,1e-10);%ԭ����Ԫ�ر�Ϊ1e-10
            U = U./repmat(norms,mFea,1);
            V = V.*repmat(norms,nSmp,1);
        end
    end

end    