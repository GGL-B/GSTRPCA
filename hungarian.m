function [C,T]=hungarian(A)
%HUNGARIAN Solve the Assignment problem using the Hungarian method.使用Hungarian 方法解决赋值问题
%
%[C,T]=hungarian(A)
%A - a square cost matrix.平方成本矩阵。
%C - the optimal assignment.最佳分配
%T - the cost of the optimal assignment.最优分配的成本。
%s.t. T = trace(A(C,:)) is minimized over all possible assignments.在所有可能的分配中最小化。

% Adapted from the FORTRAN IV code in Carpaneto and Toth, "Algorithm 548:
% Solution of the assignment problem [H]", ACM Transactions on
% Mathematical Software, 6(1):104-111, 1980.

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.
%                 Department of Computing Science, Ume? University,
%                 Sweden. 
%                 All standard disclaimers apply.

% A substantial effort was put into this code. If you use it for a
% publication or otherwise, please include an acknowledgement or at least
% notify me by email. /Niclas

[m,n]=size(A);

if (m~=n)
    error('HUNGARIAN: Cost matrix must be square!');
end

% Save original cost matrix.保存原始成本矩阵
orig=A;

% Reduce matrix.减少矩阵。
A=hminired(A);

% Do an initial assignment.完成初始任务。
[A,C,U]=hminiass(A);

% Repeat while we have unassigned rows.当我们有未分配的行时重复。
while (U(n+1))
    % Start with no path, no unchecked zeros, and no unexplored rows.开始时没有路径，没有未选中的零，也没有未探索的行。
    LR=zeros(1,n);
    LC=zeros(1,n);
    CH=zeros(1,n);
    RH=[zeros(1,n) -1];
    
    % No labelled columns.没有标签列
    SLC=[];
    
    % Start path in first unassigned row.第一个未分配行中的起始路径。
    r=U(n+1);
    % Mark row with end-of-path label.用路径结束标签标记行
    LR(r)=-1;
    % Insert row first in labelled row set.在标记的行集中首先插入行。
    SLR=r;
    
    % Repeat until we manage to find an assignable zero.重复，直到找到一个可赋值的零。
    while (1)
        % If there are free zeros in row r 如果r行中有自由零
        if (A(r,n+1)~=0)
            % ...get column of first free zero.获取第一个自由零的列。
            l=-A(r,n+1);
            
            % If there are more free zeros in row r and row r in not
            % yet marked as unexplored..如果第r行和第r行中有更多的空闲零尚未标记为未探索。。
            if (A(r,l)~=0 & RH(r)==0)
                % Insert row r first in unexplored list.在未探索列表中首先插入第r行。
                RH(r)=RH(n+1);
                RH(n+1)=r;
                
                % Mark in which column the next unexplored zero in this row
                % is. 标记此行中下一个未探索的零位于哪一列。
                CH(r)=-A(r,l);
            end
        else
            % If all rows are explored..如果浏览了所有行
            if (RH(n+1)<=0)
                % Reduce matrix.减少矩阵。
                [A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR);
            end
            
            % Re-start with first unexplored row. 从第一个未探索的行重新开始。
            r=RH(n+1);
            % Get column of next free zero in row r.获取r行中下一个空闲零的列。
            l=CH(r);
            % Advance "column of next free zero".
            CH(r)=-A(r,l);
            % If this zero is last in the list..如果此零是列表中的最后一个。
            if (A(r,l)==0)
                % ...remove row r from unexplored list.从未探索列表中删除行r。
                RH(n+1)=RH(r);
                RH(r)=0;
            end
        end
        
        % While the column l is labelled, i.e. in path.当列l被标记时，即在路径中。
        while (LC(l)~=0)
            % If row r is explored..如果探索第r行。
            if (RH(r)==0)
                % If all rows are explored..如果浏览了所有行。
                if (RH(n+1)<=0)
                    % Reduce cost matrix.降低成本矩阵。
                    [A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR);
                end
                
                % Re-start with first unexplored row.从第一个未探索的行重新开始。
                r=RH(n+1);
            end
            
            % Get column of next free zero in row r.获取r行中下一个空闲零的列
            l=CH(r);
            
            % Advance "column of next free zero".前进“下一个自由零点的列”。
            CH(r)=-A(r,l);
            
            % If this zero is last in list..如果此零是列表中的最后一个。
            if(A(r,l)==0)
                % ...remove row r from unexplored list.从未探索列表中删除行r
                RH(n+1)=RH(r);
                RH(r)=0;
            end
        end
        
        % If the column found is unassigned..如果找到的列未分配。。
        if (C(l)==0)
            % Flip all zeros along the path in LR,LC.  沿LR、LC中的路径翻转所有零。
            [A,C,U]=hmflip(A,C,LC,LR,U,l,r);
            % ...and exit to continue with next unassigned row.并退出以继续下一个未分配的行。
            break;
        else
            % ...else add zero to path.否则将零添加到路径。
            
            % Label column l with row r.用r行标记列l。
            LC(l)=r;
            
            % Add l to the set of labelled columns.将l添加到标记列集合。
            SLC=[SLC l];
            
            % Continue with the row assigned to column l.继续指定给列l的行。
            r=C(l);
            
            % Label row r with column l.用列l标记行r。
            LR(r)=l;
            
            % Add r to the set of labelled rows.将r添加到标记的行集合。
            SLR=[SLR r];
        end
    end
end

% Calculate the total cost.计算总成本。
T=sum(orig(logical(sparse(C,1:size(orig,2),1))));


function A=hminired(A)
%HMINIRED Initial reduction of cost matrix for the Hungarian method.  Hungarian方法成本矩阵的初始缩减。
%
%B=assredin(A)
%A - the unreduced cost matris.未缩减的成本矩阵。
%B - the reduced cost matrix with linked zeros in each row.在每行中具有链接零的降低成本矩阵。

% v1.0  96-06-13. Niclas Borlin, niclas@cs.umu.se.

[m,n]=size(A);

% Subtract column-minimum values from each column.从每列中减去列最小值。
colMin=min(A);
A=A-colMin(ones(n,1),:);

% Subtract row-minimum values from each row.从每行中减去行最小值。
rowMin=min(A')';
A=A-rowMin(:,ones(1,n));

% Get positions of all zeros.获取所有零的位置。
[i,j]=find(A==0);

% Extend A to give room for row zero list header column.扩展A以为零行列表标题列留出空间。
A(1,n+1)=0;
for k=1:n
    % Get all column in this row. 获取此行中的所有列
    cols=j(k==i)';
    % Insert pointers in matrix.在矩阵中插入指针
    A(k,[n+1 cols])=[-cols 0];
end


function [A,C,U]=hminiass(A)
%HMINIASS Initial assignment of the Hungarian method. Hungarian方法的初始赋值。
%
%[B,C,U]=hminiass(A)
%A - the reduced cost matrix.降低成本矩阵
%B - the reduced cost matrix, with assigned zeros removed from lists.降低成本矩阵，从列表中删除分配的零。
%C - a vector. C(J)=I means row I is assigned to column J, 矢量C（J）＝I表示行I被分配给列J，
%              i.e. there is an assigned zero in position I,J.
%U - a vector with a linked list of unassigned rows.具有未分配行的链接列表的向量。

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

[n,np1]=size(A);

% Initalize return vectors.初始化返回向量。
C=zeros(1,n);
U=zeros(1,n+1);

% Initialize last/next zero "pointers".初始化上一个/下一个零“指针”
LZ=zeros(1,n);
NZ=zeros(1,n);

for i=1:n
    % Set j to first unassigned zero in row i.将j设置为第i行中第一个未赋值的零。
	lj=n+1;
	j=-A(i,lj);

    % Repeat until we have no more zeros (j==0) or we find a zero重复，直到没有更多的零（j==0）或找到零
	% in an unassigned column (c(j)==0).在未分配的列（c（j）（0）中。
    
	while (C(j)~=0)
		% Advance lj and j in zero list.在零列表中前进lj和j。
		lj=j;
		j=-A(i,lj);
	
		% Stop if we hit end of list.如果我们到达列表末尾，请停止。
		if (j==0)
			break;
		end
	end

	if (j~=0)
		% We found a zero in an unassigned column.我们在未分配的列中找到一个零。
		
		% Assign row i to column j.将行i分配给列j。
		C(j)=i;
		
		% Remove A(i,j) from unassigned zero list.从未分配的零列表中删除A（i，j）。
		A(i,lj)=A(i,j);

		% Update next/last unassigned zero pointers.更新下一个/上一个未分配的零指针。
		NZ(i)=-A(i,j);
		LZ(i)=lj;

		% Indicate A(i,j) is an assigned zero.表示A（i，j）为指定的零。
		A(i,j)=0;
	else
		% We found no zero in an unassigned column.在未分配的列中找不到零。

		% Check all zeros in this row.选中此行中的所有零。

		lj=n+1;
		j=-A(i,lj);
		
		% Check all zeros in this row for a suitable zero in another row.检查该行中的所有零，以便在另一行中找到合适的零。
		while (j~=0)
			% Check the in the row assigned to this column.选中分配给此列的行中的。
			r=C(j);
			
			% Pick up last/next pointers.拾取上一个/下一个指针
			lm=LZ(r);
			m=NZ(r);
			
			% Check all unchecked zeros in free list of this row.选中此行自由列表中所有未选中的零。
			while (m~=0)
				% Stop if we find an unassigned column.如果找到未分配的列，请停止。
				if (C(m)==0)
					break;
				end
				
				% Advance one step in list.在列表中前进一步。
				lm=m;
				m=-A(r,lm);
			end
			
			if (m==0)
				% We failed on row r. Continue with next zero on row i.我们在第r行失败。继续第i行的下一个零。
				lj=j;
				j=-A(i,lj);
			else
				% We found a zero in an unassigned column.我们在未分配的列中找到一个零。
			
				% Replace zero at (r,m) in unassigned list with zero at (r,j) 将未分配列表中（r，m）处的零替换为（r，j）处的0
				A(r,lm)=-j;
				A(r,j)=A(r,m);
			
				% Update last/next pointers in row r.更新r行中的上一个/下一个指针。
				NZ(r)=-A(r,m);
				LZ(r)=j;
			
				% Mark A(r,m) as an assigned zero in the matrix . . .将A（r，m）标记为矩阵中指定的零。
				A(r,m)=0;
			
				% ...and in the assignment vector.以及在分配向量中。
				C(m)=r;
			
				% Remove A(i,j) from unassigned list.从未分配列表中删除A（i，j）。
				A(i,lj)=A(i,j);
			
				% Update last/next pointers in row r.更新r行中的上一个/下一个指针。
				NZ(i)=-A(i,j);
				LZ(i)=lj;
			
				% Mark A(r,m) as an assigned zero in the matrix . . .将A（r，m）标记为矩阵中指定的零。
				A(i,j)=0;
			
				% ...and in the assignment vector.以及在分配向量中。
				C(j)=i;
				
				% Stop search.停止搜索。
				break;
			end
		end
	end
end

% Create vector with list of unassigned rows.使用未分配行列表创建向量。

% Mark all rows have assignment.标记所有行都有赋值。
r=zeros(1,n);
rows=C(C~=0);
r(rows)=rows;
empty=find(r==0);

% Create vector with linked list of unassigned rows.使用未分配行的链接列表创建向量。
U=zeros(1,n+1);
U([n+1 empty])=[empty 0];


function [A,C,U]=hmflip(A,C,LC,LR,U,l,r)
%HMFLIP Flip assignment state of all zeros along a path.HMFLIP翻转路径上所有零的赋值状态。
%
%[A,C,U]=hmflip(A,C,LC,LR,U,l,r)
%Input:
%A   - the cost matrix.成本矩阵。
%C   - the assignment vector.分配向量。
%LC  - the column label vector.列标签向量。
%LR  - the row label vector.行标签向量
%U   - the 
%r,l - position of last zero in path.路径中最后一个零的位置。
%Output:
%A   - updated cost matrix.更新成本矩阵
%C   - updated assignment vector.更新的分配向量。
%U   - updated unassigned row list vector.更新了未分配的行列表矢量。

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

n=size(A,1);

while (1)
    % Move assignment in column l to row r.将l列中的赋值移动到r行
    C(l)=r;
    
    % Find zero to be removed from zero list..查找要从零列表中删除的零。
    
    % Find zero before this.在此之前查找零。
    m=find(A(r,:)==-l);
    
    % Link past this zero.链接超过此零点。
    A(r,m)=A(r,l);
    
    A(r,l)=0;
    
    % If this was the first zero of the path..如果这是路径的第一个零。
    if (LR(r)<0)
        ...remove row from unassigned row list and return. 从未分配行列表中删除行并返回。
        U(n+1)=U(r);
        U(r)=0;
        return;
    else
        
        % Move back in this row along the path and get column of next zero.
        % 沿路径在此行中向后移动并获取下一个零的列。
        l=LR(r);
        
        % Insert zero at (r,l) first in zero list.在零列表的第一个位置（r，l）插入零。
        A(r,l)=A(r,n+1);
        A(r,n+1)=-l;
        
        % Continue back along the column to get row of next zero in path.继续沿列返回，以获取路径中下一个零的行。
        r=LC(l);
    end
end


function [A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR)
%HMREDUCE Reduce parts of cost matrix in the Hungerian method. 用 Hungerian 方法减少部分成本矩阵。
%
%[A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR)
%Input:
%A   - Cost matrix.成本矩阵。
%CH  - vector of column of 'next zeros' in each row.每行中“下一个零”列的向量。
%RH  - vector with list of unexplored rows.带有未探索行列表的矢量。
%LC  - column labels.列标签。
%RC  - row labels.行标签。
%SLC - set of column labels.一组列标签。
%SLR - set of row labels.行标签集。
%
%Output:
%A   - Reduced cost matrix.降低成本矩阵。
%CH  - Updated vector of 'next zeros' in each row.更新了每行中“下一个零”的矢量
%RH  - Updated vector of unexplored rows.更新了未探索行的矢量。

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

n=size(A,1);

% Find which rows are covered, i.e. unlabelled.查找覆盖的行
coveredRows=LR==0;

% Find which columns are covered, i.e. labelled.查找覆盖的列
coveredCols=LC~=0;

r=find(~coveredRows);
c=find(~coveredCols);

% Get minimum of uncovered elements.获取最小的未覆盖元素。
m=min(min(A(r,c)));

% Subtract minimum from all uncovered elements.从所有未覆盖元素中减去最小值。
A(r,c)=A(r,c)-m;

% Check all uncovered columns..检查所有未覆盖的列。
for j=c
    % ...and uncovered rows in path order..以及按路径顺序未覆盖的行。
    for i=SLR
        % If this is a (new) zero..如果这是一个（新的）零。
        if (A(i,j)==0)
            % If the row is not in unexplored list..如果该行不在未探索列表中。
            if (RH(i)==0)
                % ...insert it first in unexplored list.在未探索列表中首先插入它。
                RH(i)=RH(n+1);
                RH(n+1)=i;
                % Mark this zero as "next free" in this row.将此零标记为此行中的“下一个空闲”。
                CH(i)=j;
            end
            % Find last unassigned zero on row I.查找行I上最后一个未分配的零。
            row=A(i,:);
            colsInList=-row(row<0);
            if (length(colsInList)==0)
                % No zeros in the list.列表中没有零。
                l=n+1;
            else
                l=colsInList(row(colsInList)==0);
            end
            % Append this zero to end of list.将此零追加到列表末尾。
            A(i,l)=-j;
        end
    end
end

% Add minimum to all doubly covered elements.将最小值添加到所有双覆盖元素。
r=find(coveredRows);
c=find(coveredCols);

% Take care of the zeros we will remove.注意我们要删除的零。
[i,j]=find(A(r,c)<=0);

i=r(i);
j=c(j);

for k=1:length(i)
    % Find zero before this in this row.在此行中，在此之前查找零。
    lj=find(A(i(k),:)==-j(k));
    % Link past it.
    A(i(k),lj)=A(i(k),j(k));
    % Mark it as assigned.将其标记为已分配
    A(i(k),j(k))=0;
end

A(r,c)=A(r,c)+m;
