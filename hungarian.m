function [C,T]=hungarian(A)
%HUNGARIAN Solve the Assignment problem using the Hungarian method.ʹ��Hungarian ���������ֵ����
%
%[C,T]=hungarian(A)
%A - a square cost matrix.ƽ���ɱ�����
%C - the optimal assignment.��ѷ���
%T - the cost of the optimal assignment.���ŷ���ĳɱ���
%s.t. T = trace(A(C,:)) is minimized over all possible assignments.�����п��ܵķ�������С����

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

% Save original cost matrix.����ԭʼ�ɱ�����
orig=A;

% Reduce matrix.���پ���
A=hminired(A);

% Do an initial assignment.��ɳ�ʼ����
[A,C,U]=hminiass(A);

% Repeat while we have unassigned rows.��������δ�������ʱ�ظ���
while (U(n+1))
    % Start with no path, no unchecked zeros, and no unexplored rows.��ʼʱû��·����û��δѡ�е��㣬Ҳû��δ̽�����С�
    LR=zeros(1,n);
    LC=zeros(1,n);
    CH=zeros(1,n);
    RH=[zeros(1,n) -1];
    
    % No labelled columns.û�б�ǩ��
    SLC=[];
    
    % Start path in first unassigned row.��һ��δ�������е���ʼ·����
    r=U(n+1);
    % Mark row with end-of-path label.��·��������ǩ�����
    LR(r)=-1;
    % Insert row first in labelled row set.�ڱ�ǵ��м������Ȳ����С�
    SLR=r;
    
    % Repeat until we manage to find an assignable zero.�ظ���ֱ���ҵ�һ���ɸ�ֵ���㡣
    while (1)
        % If there are free zeros in row r ���r������������
        if (A(r,n+1)~=0)
            % ...get column of first free zero.��ȡ��һ����������С�
            l=-A(r,n+1);
            
            % If there are more free zeros in row r and row r in not
            % yet marked as unexplored..�����r�к͵�r�����и���Ŀ�������δ���Ϊδ̽������
            if (A(r,l)~=0 & RH(r)==0)
                % Insert row r first in unexplored list.��δ̽���б������Ȳ����r�С�
                RH(r)=RH(n+1);
                RH(n+1)=r;
                
                % Mark in which column the next unexplored zero in this row
                % is. ��Ǵ�������һ��δ̽������λ����һ�С�
                CH(r)=-A(r,l);
            end
        else
            % If all rows are explored..��������������
            if (RH(n+1)<=0)
                % Reduce matrix.���پ���
                [A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR);
            end
            
            % Re-start with first unexplored row. �ӵ�һ��δ̽���������¿�ʼ��
            r=RH(n+1);
            % Get column of next free zero in row r.��ȡr������һ����������С�
            l=CH(r);
            % Advance "column of next free zero".
            CH(r)=-A(r,l);
            % If this zero is last in the list..����������б��е����һ����
            if (A(r,l)==0)
                % ...remove row r from unexplored list.��δ̽���б���ɾ����r��
                RH(n+1)=RH(r);
                RH(r)=0;
            end
        end
        
        % While the column l is labelled, i.e. in path.����l�����ʱ������·���С�
        while (LC(l)~=0)
            % If row r is explored..���̽����r�С�
            if (RH(r)==0)
                % If all rows are explored..�������������С�
                if (RH(n+1)<=0)
                    % Reduce cost matrix.���ͳɱ�����
                    [A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR);
                end
                
                % Re-start with first unexplored row.�ӵ�һ��δ̽���������¿�ʼ��
                r=RH(n+1);
            end
            
            % Get column of next free zero in row r.��ȡr������һ�����������
            l=CH(r);
            
            % Advance "column of next free zero".ǰ������һ�����������С���
            CH(r)=-A(r,l);
            
            % If this zero is last in list..����������б��е����һ����
            if(A(r,l)==0)
                % ...remove row r from unexplored list.��δ̽���б���ɾ����r
                RH(n+1)=RH(r);
                RH(r)=0;
            end
        end
        
        % If the column found is unassigned..����ҵ�����δ���䡣��
        if (C(l)==0)
            % Flip all zeros along the path in LR,LC.  ��LR��LC�е�·����ת�����㡣
            [A,C,U]=hmflip(A,C,LC,LR,U,l,r);
            % ...and exit to continue with next unassigned row.���˳��Լ�����һ��δ������С�
            break;
        else
            % ...else add zero to path.��������ӵ�·����
            
            % Label column l with row r.��r�б����l��
            LC(l)=r;
            
            % Add l to the set of labelled columns.��l��ӵ�����м��ϡ�
            SLC=[SLC l];
            
            % Continue with the row assigned to column l.����ָ������l���С�
            r=C(l);
            
            % Label row r with column l.����l�����r��
            LR(r)=l;
            
            % Add r to the set of labelled rows.��r��ӵ���ǵ��м��ϡ�
            SLR=[SLR r];
        end
    end
end

% Calculate the total cost.�����ܳɱ���
T=sum(orig(logical(sparse(C,1:size(orig,2),1))));


function A=hminired(A)
%HMINIRED Initial reduction of cost matrix for the Hungarian method.  Hungarian�����ɱ�����ĳ�ʼ������
%
%B=assredin(A)
%A - the unreduced cost matris.δ�����ĳɱ�����
%B - the reduced cost matrix with linked zeros in each row.��ÿ���о���������Ľ��ͳɱ�����

% v1.0  96-06-13. Niclas Borlin, niclas@cs.umu.se.

[m,n]=size(A);

% Subtract column-minimum values from each column.��ÿ���м�ȥ����Сֵ��
colMin=min(A);
A=A-colMin(ones(n,1),:);

% Subtract row-minimum values from each row.��ÿ���м�ȥ����Сֵ��
rowMin=min(A')';
A=A-rowMin(:,ones(1,n));

% Get positions of all zeros.��ȡ�������λ�á�
[i,j]=find(A==0);

% Extend A to give room for row zero list header column.��չA��Ϊ�����б�����������ռ䡣
A(1,n+1)=0;
for k=1:n
    % Get all column in this row. ��ȡ�����е�������
    cols=j(k==i)';
    % Insert pointers in matrix.�ھ����в���ָ��
    A(k,[n+1 cols])=[-cols 0];
end


function [A,C,U]=hminiass(A)
%HMINIASS Initial assignment of the Hungarian method. Hungarian�����ĳ�ʼ��ֵ��
%
%[B,C,U]=hminiass(A)
%A - the reduced cost matrix.���ͳɱ�����
%B - the reduced cost matrix, with assigned zeros removed from lists.���ͳɱ����󣬴��б���ɾ��������㡣
%C - a vector. C(J)=I means row I is assigned to column J, ʸ��C��J����I��ʾ��I���������J��
%              i.e. there is an assigned zero in position I,J.
%U - a vector with a linked list of unassigned rows.����δ�����е������б��������

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

[n,np1]=size(A);

% Initalize return vectors.��ʼ������������
C=zeros(1,n);
U=zeros(1,n+1);

% Initialize last/next zero "pointers".��ʼ����һ��/��һ���㡰ָ�롱
LZ=zeros(1,n);
NZ=zeros(1,n);

for i=1:n
    % Set j to first unassigned zero in row i.��j����Ϊ��i���е�һ��δ��ֵ���㡣
	lj=n+1;
	j=-A(i,lj);

    % Repeat until we have no more zeros (j==0) or we find a zero�ظ���ֱ��û�и�����㣨j==0�����ҵ���
	% in an unassigned column (c(j)==0).��δ������У�c��j����0���С�
    
	while (C(j)~=0)
		% Advance lj and j in zero list.�����б���ǰ��lj��j��
		lj=j;
		j=-A(i,lj);
	
		% Stop if we hit end of list.������ǵ����б�ĩβ����ֹͣ��
		if (j==0)
			break;
		end
	end

	if (j~=0)
		% We found a zero in an unassigned column.������δ����������ҵ�һ���㡣
		
		% Assign row i to column j.����i�������j��
		C(j)=i;
		
		% Remove A(i,j) from unassigned zero list.��δ��������б���ɾ��A��i��j����
		A(i,lj)=A(i,j);

		% Update next/last unassigned zero pointers.������һ��/��һ��δ�������ָ�롣
		NZ(i)=-A(i,j);
		LZ(i)=lj;

		% Indicate A(i,j) is an assigned zero.��ʾA��i��j��Ϊָ�����㡣
		A(i,j)=0;
	else
		% We found no zero in an unassigned column.��δ����������Ҳ����㡣

		% Check all zeros in this row.ѡ�д����е������㡣

		lj=n+1;
		j=-A(i,lj);
		
		% Check all zeros in this row for a suitable zero in another row.�������е������㣬�Ա�����һ�����ҵ����ʵ��㡣
		while (j~=0)
			% Check the in the row assigned to this column.ѡ�з�������е����еġ�
			r=C(j);
			
			% Pick up last/next pointers.ʰȡ��һ��/��һ��ָ��
			lm=LZ(r);
			m=NZ(r);
			
			% Check all unchecked zeros in free list of this row.ѡ�д��������б�������δѡ�е��㡣
			while (m~=0)
				% Stop if we find an unassigned column.����ҵ�δ������У���ֹͣ��
				if (C(m)==0)
					break;
				end
				
				% Advance one step in list.���б���ǰ��һ����
				lm=m;
				m=-A(r,lm);
			end
			
			if (m==0)
				% We failed on row r. Continue with next zero on row i.�����ڵ�r��ʧ�ܡ�������i�е���һ���㡣
				lj=j;
				j=-A(i,lj);
			else
				% We found a zero in an unassigned column.������δ����������ҵ�һ���㡣
			
				% Replace zero at (r,m) in unassigned list with zero at (r,j) ��δ�����б��У�r��m���������滻Ϊ��r��j������0
				A(r,lm)=-j;
				A(r,j)=A(r,m);
			
				% Update last/next pointers in row r.����r���е���һ��/��һ��ָ�롣
				NZ(r)=-A(r,m);
				LZ(r)=j;
			
				% Mark A(r,m) as an assigned zero in the matrix . . .��A��r��m�����Ϊ������ָ�����㡣
				A(r,m)=0;
			
				% ...and in the assignment vector.�Լ��ڷ��������С�
				C(m)=r;
			
				% Remove A(i,j) from unassigned list.��δ�����б���ɾ��A��i��j����
				A(i,lj)=A(i,j);
			
				% Update last/next pointers in row r.����r���е���һ��/��һ��ָ�롣
				NZ(i)=-A(i,j);
				LZ(i)=lj;
			
				% Mark A(r,m) as an assigned zero in the matrix . . .��A��r��m�����Ϊ������ָ�����㡣
				A(i,j)=0;
			
				% ...and in the assignment vector.�Լ��ڷ��������С�
				C(j)=i;
				
				% Stop search.ֹͣ������
				break;
			end
		end
	end
end

% Create vector with list of unassigned rows.ʹ��δ�������б���������

% Mark all rows have assignment.��������ж��и�ֵ��
r=zeros(1,n);
rows=C(C~=0);
r(rows)=rows;
empty=find(r==0);

% Create vector with linked list of unassigned rows.ʹ��δ�����е������б���������
U=zeros(1,n+1);
U([n+1 empty])=[empty 0];


function [A,C,U]=hmflip(A,C,LC,LR,U,l,r)
%HMFLIP Flip assignment state of all zeros along a path.HMFLIP��ת·����������ĸ�ֵ״̬��
%
%[A,C,U]=hmflip(A,C,LC,LR,U,l,r)
%Input:
%A   - the cost matrix.�ɱ�����
%C   - the assignment vector.����������
%LC  - the column label vector.�б�ǩ������
%LR  - the row label vector.�б�ǩ����
%U   - the 
%r,l - position of last zero in path.·�������һ�����λ�á�
%Output:
%A   - updated cost matrix.���³ɱ�����
%C   - updated assignment vector.���µķ���������
%U   - updated unassigned row list vector.������δ��������б�ʸ����

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

n=size(A,1);

while (1)
    % Move assignment in column l to row r.��l���еĸ�ֵ�ƶ���r��
    C(l)=r;
    
    % Find zero to be removed from zero list..����Ҫ�����б���ɾ�����㡣
    
    % Find zero before this.�ڴ�֮ǰ�����㡣
    m=find(A(r,:)==-l);
    
    % Link past this zero.���ӳ�������㡣
    A(r,m)=A(r,l);
    
    A(r,l)=0;
    
    % If this was the first zero of the path..�������·���ĵ�һ���㡣
    if (LR(r)<0)
        ...remove row from unassigned row list and return. ��δ�������б���ɾ���в����ء�
        U(n+1)=U(r);
        U(r)=0;
        return;
    else
        
        % Move back in this row along the path and get column of next zero.
        % ��·���ڴ���������ƶ�����ȡ��һ������С�
        l=LR(r);
        
        % Insert zero at (r,l) first in zero list.�����б�ĵ�һ��λ�ã�r��l�������㡣
        A(r,l)=A(r,n+1);
        A(r,n+1)=-l;
        
        % Continue back along the column to get row of next zero in path.�������з��أ��Ի�ȡ·������һ������С�
        r=LC(l);
    end
end


function [A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR)
%HMREDUCE Reduce parts of cost matrix in the Hungerian method. �� Hungerian �������ٲ��ֳɱ�����
%
%[A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR)
%Input:
%A   - Cost matrix.�ɱ�����
%CH  - vector of column of 'next zeros' in each row.ÿ���С���һ���㡱�е�������
%RH  - vector with list of unexplored rows.����δ̽�����б��ʸ����
%LC  - column labels.�б�ǩ��
%RC  - row labels.�б�ǩ��
%SLC - set of column labels.һ���б�ǩ��
%SLR - set of row labels.�б�ǩ����
%
%Output:
%A   - Reduced cost matrix.���ͳɱ�����
%CH  - Updated vector of 'next zeros' in each row.������ÿ���С���һ���㡱��ʸ��
%RH  - Updated vector of unexplored rows.������δ̽���е�ʸ����

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

n=size(A,1);

% Find which rows are covered, i.e. unlabelled.���Ҹ��ǵ���
coveredRows=LR==0;

% Find which columns are covered, i.e. labelled.���Ҹ��ǵ���
coveredCols=LC~=0;

r=find(~coveredRows);
c=find(~coveredCols);

% Get minimum of uncovered elements.��ȡ��С��δ����Ԫ�ء�
m=min(min(A(r,c)));

% Subtract minimum from all uncovered elements.������δ����Ԫ���м�ȥ��Сֵ��
A(r,c)=A(r,c)-m;

% Check all uncovered columns..�������δ���ǵ��С�
for j=c
    % ...and uncovered rows in path order..�Լ���·��˳��δ���ǵ��С�
    for i=SLR
        % If this is a (new) zero..�������һ�����µģ��㡣
        if (A(i,j)==0)
            % If the row is not in unexplored list..������в���δ̽���б��С�
            if (RH(i)==0)
                % ...insert it first in unexplored list.��δ̽���б������Ȳ�������
                RH(i)=RH(n+1);
                RH(n+1)=i;
                % Mark this zero as "next free" in this row.��������Ϊ�����еġ���һ�����С���
                CH(i)=j;
            end
            % Find last unassigned zero on row I.������I�����һ��δ������㡣
            row=A(i,:);
            colsInList=-row(row<0);
            if (length(colsInList)==0)
                % No zeros in the list.�б���û���㡣
                l=n+1;
            else
                l=colsInList(row(colsInList)==0);
            end
            % Append this zero to end of list.������׷�ӵ��б�ĩβ��
            A(i,l)=-j;
        end
    end
end

% Add minimum to all doubly covered elements.����Сֵ��ӵ�����˫����Ԫ�ء�
r=find(coveredRows);
c=find(coveredCols);

% Take care of the zeros we will remove.ע������Ҫɾ�����㡣
[i,j]=find(A(r,c)<=0);

i=r(i);
j=c(j);

for k=1:length(i)
    % Find zero before this in this row.�ڴ����У��ڴ�֮ǰ�����㡣
    lj=find(A(i(k),:)==-j(k));
    % Link past it.
    A(i(k),lj)=A(i(k),j(k));
    % Mark it as assigned.������Ϊ�ѷ���
    A(i(k),j(k))=0;
end

A(r,c)=A(r,c)+m;
