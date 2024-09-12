function [W, elapse] = constructW(fea,options)
%	Usage:
%	W = constructW(fea,options)
%
%	fea: Rows of vectors of data points. Each row is x_i 数据点矢量行。每行为x_i

%	%每行对应一个点,列对应不同属性值

%   options: Struct value in Matlab. The fields in options that can be set:Matlab中的结构值。选项中可设置的字段：
%           Metric -  Choices are:
%               'Euclidean' - Will use the Euclidean distance of two data 
%                             points to evaluate the "closeness" between 
%                             them. [Default One]我们将使用两个数据点的欧几里德距离来评估两者之间的“接近度”
%               'Cosine'    - Will use the cosine value of two vectors
%                             to evaluate the "closeness" between them.
%                             A popular similarity measure used in
%                             Information Retrieval.我将使用两个向量的余弦值来评估它们之间的“接近度”。使用一种流行的相似性度量
%                  
%           NeighborMode -  Indicates how to construct the graph. Choices指示如何构造图形
%                           are: [Default 'KNN']
%                'KNN'            -  k = 0
%                                       Complete graph
%                                    k > 0
%                                      Put an edge between two nodes if and
%                                      only if they are among the k nearst
%                                      neighbors of each other. You are
%                                      required to provide the parameter k in
%                                      the options. Default k=5. 在两个节点之间放置一条边，当且仅当它们位于彼此的k个最近邻居之间。需要在选项中提供参数k。默认值k=5
%               'Supervised'      -  k = 0
%                                       Put an edge between two nodes if and
%                                       only if they belong to same class. 当且仅当两个节点属于同一类时，在两个节点之间切出一条边。
%                                    k > 0
%                                       Put an edge between two nodes if
%                                       they belong to same class and they
%                                       are among the k nearst neighbors of
%                                       each other. 如果两个节点属于同一类，并且它们在彼此的k个最近邻居之间，则在两个节点之间创建一条边。
%                                    Default: k=0
%                                   You are required to provide the label
%                                   information gnd in the options.您需要提供标签选项中的信息gnd。
%                                              
%           WeightMode   -  Indicates how to assign weights for each edge
%                           in the graph. Choices are:   指示如何为图形中的每个分配权重
%               'Binary'       - 0-1 weighting. Every edge receiveds weight
%                                of 1. [Default One]  每条边接收1的权重
%               'HeatKernel'   - If nodes i and j are connected, put weight
%                                W_ij = exp(-norm(x_i - x_j)/2t^2). This
%                                weight mode can only be used under
%                                'Euclidean' metric and you are required to
%                                provide the parameter t. 如果节点i和j连接，则将权重W_ij=exp（-范数（x_i-x_j）/2t^2）。此权重模式只能在“欧几里得”度量下使用，
%               'Cosine'       - If nodes i and j are connected, put weight
%                                cosine(x_i,x_j). Can only be used under
%                                'Cosine' metric. 如果节点i和j是连接的，则输入权重余弦（x_i，x_j）。只能在“余弦”度量下使用。
%               
%            k         -   The parameter needed under 'KNN' NeighborMode.
%                          Default will be 5.“KNN”模式默认值下所需的参数为5。
%            gnd       -   The parameter needed under 'Supervised'
%                          NeighborMode.  Colunm vector of the label
%                          information for each data point.受监督”邻居模式下需要的参数。每个数据点的标签信息的列向量。
%            bLDA      -   0 or 1. Only effective under 'Supervised'
%                          NeighborMode. If 1, the graph will be constructed
%                          to make LPP exactly same as LDA. Default will be
%                          0. 
%            t         -   The parameter needed under 'HeatKernel'
%                          WeightMode. Default will be 1
%         bNormalized  -   0 or 1. Only effective under 'Cosine' metric.
%                          Indicates whether the fea are already be
%                          normalized to 1. Default will be 0
%      bSelfConnected  -   0 or 1. Indicates whether W(i,i) == 1. Default 1
%                          if 'Supervised' NeighborMode & bLDA == 1,
%                          bSelfConnected will always be 1. Default 1.
%“HeatKernel”权重模式下需要的参数。默认值为0或1,
%
%    Examples:
%
%       fea = rand(50,15);
%       options = [];
%       options.Metric = 'Euclidean';
%       options.NeighborMode = 'KNN';
%       options.k = 5;
%       options.WeightMode = 'HeatKernel';
%       options.t = 1;
%       W = constructW(fea,options);
%       
%       
%       fea = rand(50,15);
%       gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];
%       options = [];
%       options.Metric = 'Euclidean';
%       options.NeighborMode = 'Supervised';
%       options.gnd = gnd;
%       options.WeightMode = 'HeatKernel';
%       options.t = 1;
%       W = constructW(fea,options);
%       
%       
%       fea = rand(50,15);
%       gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];
%       options = [];
%       options.Metric = 'Euclidean';
%       options.NeighborMode = 'Supervised';
%       options.gnd = gnd;
%       options.bLDA = 1;
%       W = constructW(fea,options);      
%       
%
%    For more details about the different ways to construct the W, please
%    refer:
%       Deng Cai, Xiaofei He and Jiawei Han, "Document Clustering Using
%       Locality Preserving Indexing" IEEE TKDE, Dec. 2005.
%    
%
%    Written by Deng Cai (dengcai2 AT cs.uiuc.edu), April/2004, Feb/2006,
%                                             May/2007
% 

if (~exist('options','var'))
   options = [];
else
   if ~isstruct(options) 
       error('parameter error!');
   end
end

%=================================================
if ~isfield(options,'Metric')
    options.Metric = 'Cosine';
end

switch lower(options.Metric)
    case {lower('Euclidean')}
    case {lower('Cosine')}
        if ~isfield(options,'bNormalized')
            options.bNormalized = 0;
        end
    otherwise
        error('Metric does not exist!');
end

%=================================================
if ~isfield(options,'NeighborMode')
    options.NeighborMode = 'KNN';
end

switch lower(options.NeighborMode)
    case {lower('KNN')}  %For simplicity, we include the data point itself in the kNN  为了简单起见，我们将数据点本身包含在kNN中
        if ~isfield(options,'k')
            options.k = 5;
        end
    case {lower('Supervised')}
        if ~isfield(options,'bLDA')
            options.bLDA = 0;
        end
        if options.bLDA
            options.bSelfConnected = 1;
        end
        if ~isfield(options,'k')
            options.k = 0;
        end
        if ~isfield(options,'gnd')
            error('Label(gnd) should be provided under ''Supervised'' NeighborMode!');
        end
        if ~isempty(fea) && length(options.gnd) ~= size(fea,1)
            error('gnd doesn''t match with fea!');
        end
    otherwise
        error('NeighborMode does not exist!');
end

%=================================================

if ~isfield(options,'WeightMode')
    options.WeightMode = 'Binary';
end

bBinary = 0;
switch lower(options.WeightMode)
    case {lower('Binary')}
        bBinary = 1; 
    case {lower('HeatKernel')}
        if ~strcmpi(options.Metric,'Euclidean')
            warning('''HeatKernel'' WeightMode should be used under ''Euclidean'' Metric!');
            options.Metric = 'Euclidean';
        end
        if ~isfield(options,'t')
            options.t = 1;
        end
    case {lower('Cosine')}
        if ~strcmpi(options.Metric,'Cosine')
            warning('''Cosine'' WeightMode should be used under ''Cosine'' Metric!');
            options.Metric = 'Cosine';
        end
        if ~isfield(options,'bNormalized')
            options.bNormalized = 0;
        end
    otherwise
        error('WeightMode does not exist!');
end

%=================================================

if ~isfield(options,'bSelfConnected')
    options.bSelfConnected = 1;
end

%=================================================
tmp_T = cputime;

if isfield(options,'gnd') 
    nSmp = length(options.gnd);
else
    nSmp = size(fea,1);
end
maxM = 62500000; %500M
BlockSize = floor(maxM/(nSmp*3));


if strcmpi(options.NeighborMode,'Supervised')
    Label = unique(options.gnd);
    nLabel = length(Label);
    if options.bLDA
        G = zeros(nSmp,nSmp);
        for idx=1:nLabel
            classIdx = options.gnd==Label(idx);
            G(classIdx,classIdx) = 1/sum(classIdx);
        end
        W = sparse(G);
        elapse = cputime - tmp_T;
        return;
    end
    
    switch lower(options.WeightMode)
        case {lower('Binary')}
            if options.k > 0
                G = zeros(nSmp*(options.k+1),3);
                idNow = 0;
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    D = EuDist2(fea(classIdx,:),[],0);
                    [dump idx] = sort(D,2); % sort each row 对每行排序
                    clear D dump;
                    idx = idx(:,1:options.k+1);
                    
                    nSmpClass = length(classIdx)*(options.k+1);
                    G(idNow+1:nSmpClass+idNow,1) = repmat(classIdx,[options.k+1,1]);
                    G(idNow+1:nSmpClass+idNow,2) = classIdx(idx(:));
                    G(idNow+1:nSmpClass+idNow,3) = 1;
                    idNow = idNow+nSmpClass;
                    clear idx
                end
                G = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
                G = max(G,G');
            else
                G = zeros(nSmp,nSmp);
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    G(classIdx,classIdx) = 1;
                end
            end
            
            if ~options.bSelfConnected
                for i=1:size(G,1)
                    G(i,i) = 0;
                end
            end
            
            W = sparse(G);
        case {lower('HeatKernel')}
            if options.k > 0
                G = zeros(nSmp*(options.k+1),3);
                idNow = 0;
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    D = EuDist2(fea(classIdx,:),[],0);
                    [dump idx] = sort(D,2); % sort each row 对每行排序
                    clear D;
                    idx = idx(:,1:options.k+1);
                    dump = dump(:,1:options.k+1);
                    dump = exp(-dump/(2*options.t^2));
                    
                    nSmpClass = length(classIdx)*(options.k+1);
                    G(idNow+1:nSmpClass+idNow,1) = repmat(classIdx,[options.k+1,1]);
                    G(idNow+1:nSmpClass+idNow,2) = classIdx(idx(:));
                    G(idNow+1:nSmpClass+idNow,3) = dump(:);
                    idNow = idNow+nSmpClass;
                    clear dump idx
                end
                G = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
            else
                G = zeros(nSmp,nSmp);
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    D = EuDist2(fea(classIdx,:),[],0);
                    D = exp(-D/(2*options.t^2));
                    G(classIdx,classIdx) = D;
                end
            end
            
            if ~options.bSelfConnected
                for i=1:size(G,1)
                    G(i,i) = 0;
                end
            end

            W = sparse(max(G,G'));
        case {lower('Cosine')}
            if ~options.bNormalized
                [nSmp, nFea] = size(fea);
                if issparse(fea)
                    fea2 = fea';
                    feaNorm = sum(fea2.^2,1).^.5;
                    for i = 1:nSmp
                        fea2(:,i) = fea2(:,i) ./ max(1e-10,feaNorm(i));
                    end
                    fea = fea2';
                    clear fea2;
                else
                    feaNorm = sum(fea.^2,2).^.5;
                    for i = 1:nSmp
                        fea(i,:) = fea(i,:) ./ max(1e-12,feaNorm(i));
                    end
                end

            end

            if options.k > 0
                G = zeros(nSmp*(options.k+1),3);
                idNow = 0;
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    D = fea(classIdx,:)*fea(classIdx,:)';
                    [dump idx] = sort(-D,2); % sort each row 对每行排序
                    clear D;
                    idx = idx(:,1:options.k+1);
                    dump = -dump(:,1:options.k+1);
                    
                    nSmpClass = length(classIdx)*(options.k+1);
                    G(idNow+1:nSmpClass+idNow,1) = repmat(classIdx,[options.k+1,1]);
                    G(idNow+1:nSmpClass+idNow,2) = classIdx(idx(:));
                    G(idNow+1:nSmpClass+idNow,3) = dump(:);
                    idNow = idNow+nSmpClass;
                    clear dump idx
                end
                G = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
            else
                G = zeros(nSmp,nSmp);
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    G(classIdx,classIdx) = fea(classIdx,:)*fea(classIdx,:)';
                end
            end

            if ~options.bSelfConnected
                for i=1:size(G,1)
                    G(i,i) = 0;
                end
            end

            W = sparse(max(G,G'));
        otherwise
            error('WeightMode does not exist!');
    end
    elapse = cputime - tmp_T;
    return;
end


if strcmpi(options.NeighborMode,'KNN') && (options.k > 0)
    if strcmpi(options.Metric,'Euclidean')
        G = zeros(nSmp*(options.k+1),3);
        for i = 1:ceil(nSmp/BlockSize)
            if i == ceil(nSmp/BlockSize)
                smpIdx = (i-1)*BlockSize+1:nSmp;
                dist = EuDist2(fea(smpIdx,:),fea,0);
                dist = full(dist);
                [dump idx] = sort(dist,2); % sort each row 对每行排序
                idx = idx(:,1:options.k+1);
                dump = dump(:,1:options.k+1);
                if ~bBinary
                    dump = exp(-dump/(2*options.t^2));
                end
                
                G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),1) = repmat(smpIdx',[options.k+1,1]);
                G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),2) = idx(:);
                if ~bBinary
                    G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),3) = dump(:);
                else
                    G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),3) = 1;
                end
            else
                smpIdx = (i-1)*BlockSize+1:i*BlockSize;
                dist = EuDist2(fea(smpIdx,:),fea,0);
                dist = full(dist);
                [dump idx] = sort(dist,2); % sort each row 对每行排序
                idx = idx(:,1:options.k+1);
                dump = dump(:,1:options.k+1);
                if ~bBinary
                    dump = exp(-dump/(2*options.t^2));
                end
                
                G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),1) = repmat(smpIdx',[options.k+1,1]);
                G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),2) = idx(:);
                if ~bBinary
                    G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),3) = dump(:);
                else
                    G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),3) = 1;
                end
            end
        end

        W = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
    else
        if ~options.bNormalized
            [nSmp, nFea] = size(fea);
            if issparse(fea)
                fea2 = fea';
                clear fea;
                for i = 1:nSmp
                    fea2(:,i) = fea2(:,i) ./ max(1e-10,sum(fea2(:,i).^2,1).^.5);
                end
                fea = fea2';
                clear fea2;
            else
                feaNorm = sum(fea.^2,2).^.5;
                for i = 1:nSmp
                    fea(i,:) = fea(i,:) ./ max(1e-12,feaNorm(i));
                end
            end
        end
        
        G = zeros(nSmp*(options.k+1),3);
        for i = 1:ceil(nSmp/BlockSize)
            if i == ceil(nSmp/BlockSize)
                smpIdx = (i-1)*BlockSize+1:nSmp;
                dist = fea(smpIdx,:)*fea';
                dist = full(dist);
                [dump idx] = sort(-dist,2); % sort each row 对每行排序
                idx = idx(:,1:options.k+1);
                dump = -dump(:,1:options.k+1);

                G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),1) = repmat(smpIdx',[options.k+1,1]);
                G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),2) = idx(:);
                G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),3) = dump(:);
            else
                smpIdx = (i-1)*BlockSize+1:i*BlockSize;
                dist = fea(smpIdx,:)*fea';
                dist = full(dist);
                [dump idx] = sort(-dist,2); % sort each row 对每行排序
                idx = idx(:,1:options.k+1);
                dump = -dump(:,1:options.k+1);

                G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),1) = repmat(smpIdx',[options.k+1,1]);
                G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),2) = idx(:);
                G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),3) = dump(:);
            end
        end

        W = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
    end
    
    if strcmpi(options.WeightMode,'Binary')
        W(find(W)) = 1;
    end
    
    if isfield(options,'bSemiSupervised') && options.bSemiSupervised
        tmpgnd = options.gnd(options.semiSplit);
        
        Label = unique(tmpgnd);
        nLabel = length(Label);
        G = zeros(sum(options.semiSplit),sum(options.semiSplit));
        for idx=1:nLabel
            classIdx = tmpgnd==Label(idx);
            G(classIdx,classIdx) = 1;
        end
        Wsup = sparse(G);
        if ~isfield(options,'SameCategoryWeight')
            options.SameCategoryWeight = 1;
        end
        W(options.semiSplit,options.semiSplit) = (Wsup>0)*options.SameCategoryWeight;
    end
    
    if ~options.bSelfConnected
        for i=1:size(W,1)
            W(i,i) = 0;
        end
    end

    W = max(W,W');
    
    elapse = cputime - tmp_T;
    return;
end


% strcmpi(options.NeighborMode,'KNN') & (options.k == 0)
% Complete Graph 完整图表

if strcmpi(options.Metric,'Euclidean')
    W = EuDist2(fea,[],0);
    W = exp(-W/(2*options.t^2));
else
    if ~options.bNormalized
%         feaNorm = sum(fea.^2,2).^.5;
%         fea = fea ./ repmat(max(1e-10,feaNorm),1,size(fea,2));
        [nSmp, nFea] = size(fea);
        if issparse(fea)
            fea2 = fea';
            feaNorm = sum(fea2.^2,1).^.5;
            for i = 1:nSmp
                fea2(:,i) = fea2(:,i) ./ max(1e-10,feaNorm(i));
            end
            fea = fea2';
            clear fea2;
        else
            feaNorm = sum(fea.^2,2).^.5;
            for i = 1:nSmp
                fea(i,:) = fea(i,:) ./ max(1e-12,feaNorm(i));
            end
        end
    end
    
    W = full(fea*fea');
end

if ~options.bSelfConnected
    for i=1:size(W,1)
        W(i,i) = 0;
    end
end

W = max(W,W');



elapse = cputime - tmp_T;


