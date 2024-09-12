   
    function [ari,ami,nmi,accuracy] = SS(Y3)  
   
    X=Y3;
    %load('10X_inhouse')%Loading data
     load('Specter')%Loading data
%main 
    label=labels;
%     label=true_labs;
    k2=length(unique(label));% unique(label)£¬The returned value is the same as the label, but without duplicate elements. The generated result vectors are sorted in ascending order.
    scala = 0.2;  % 
    X_1 = [];
    X_2 = [];
    y_1 = [];
    y_2 = [];
    c = [];
    for labell=1:length(unique(label))
        cate = find(label==labell);%Record the position in the label that matches the label value
        half = ceil(length(cate)*scala);% scala = 0.2
%         half = int32(length(cate)*scala);
        local = randperm(length(cate),half);%Return a row vector of half elements, where the elements of this row vector are integers from 1 to length (cate)
        local_lab = cate(local);  %%Select the data at the local position from the cate
        local_non = setdiff(cate,local_lab);   % Return elements that are in the vector 'cate' but not in the vector 'local_lab'
        X_1 = [X_1,X(:,local_lab)];
        X_2 = [X_2,X(:,local_non)];
        y_1 = [y_1;label(local_lab,:)]; %%Select data from the label at the position of the local_lab line
        y_2 = [y_2;label(local_non,:)];%%Select the data at the local_non row position from the label
        c = [c;label(local_lab)];%%Select data from the label at the local_1ab position, test set: number of sample categories
    end
    X = [X_1,X_2];
    label = [y_1;y_2];
    A = labelA(label,c,k2);
    %============================
    options = [];
    option.Metric = 'Cosine';
    options.NeighborMode = 'KNN';%KNN£¬
    options.k =5;%5 nearest neighbors£¬
    options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 

    W = constructW(X',options);

    clear options;
    options = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num = length(c);%
    lambda=2;%
    k1=500; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [U_final,Z_final,B_final, F_final, S_final] = SSNMDI_model(X,A,lambda,k1,k2,W,options,num);
 
    %%%%%%%%%%% Clustering cell type label£¬
        for e=1:size(F_final,2) %%%F_final 
        v=F_final(:,e);
        ma=max(v);%ma=0.1685
        [s,t]=find(v==ma);%v=ma, position
        l(e)=s;%l(220)=6
        end
        %%%%%%%%%%%%%%Performance evaluation
        label(1:num,:)=[];%num = length(c)
        ll=label;%%%  the label originally identified by the authors¡£
        l=l';
        l(1:num,:)=[];%num = length(c),
        [newl] = bestMap(ll,l);%% Permute label of l to match ll as good as possible  
        nmi=compute_NMI(ll,newl) %% Calculating the Normalized Mutual Information (NMI)
        ami=AMI(ll,newl)%% Calculating the Adjusted Mutual Information (AMI)£¬
        ari = ARI(ll,max(ll),newl,max(newl)) %% Calculating the Adjusted Rand Index (ARI)£¬
        pre_label =newl;
        if ~isempty(ll) 
        exact = find(pre_label == ll);%ll=label,
        accuracy = length(exact)/length(newl) %% Calculating the accuracy
        else
        accuracy = []
        end
    end
