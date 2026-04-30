function P = modular_GaussElimination(A, b, q)

P= mod([A, b], q); % constructing the new augmented matrix P 
[row, col] = size( P); % Calculating the size of augmented matrix, P   
for i = 1:row-1 % Finding zeros of lower triangular matrix.
    if P(i,i) == 0 
        after_Pi = P(i+1:end, i);
        Index = find(after_Pi~=0);
        if isempty(Index)==1
            continue;
        end
        temp = P(i, :);
        
        P(i, :) = P(i+Index(1), :);
        P(i+Index(1), :) = temp;
    end
    a=P(i,i);  
    P(i,:)= mod(P(i,:)*minv(a, q), q);    %P(i,:)= P(i,:)/a;
    for j=i+1:row     
        P(j,:)= mod(P(j,:)- P(j,i)* P(i,:), q);
    end
end

a=P(row,row);  
if a~=0
    P(row,:)= mod(P(row,:)*minv(a, q), q); %P(row,:)= P(row,:)/a;
end
for i=row:-1:2   % Finding zeros of the upper triangular matrix.
    for j=i-1:-1:1    
        P(j,:)= mod(P(j,:)- P(j,i)* P(i,:), q);
    end
end 





 %
function pp = gauss_jordan(A, b)
% %Gauss elimination method [m,n]=size(A);
% clc
A= [ 2 3 6; 1 7 3; 5 1 3]; % Write the coefficient matrix, A. where  the system: AX=B.
b = [3; 9; 7]; % Write the constants matrix, B

P= [A, b]; % constructing the new augmented matrix P 
[row, col] = size( P); % Calculating the size of augmented matrix, P   
for i = 1:row-1 % Finding zeros of lower triangular matrix.
         if P(i,i) == 0 
             diag_P = diag(P);
             Index = find(diag_P(i+1:end)~=0);
             if isempty(Index)==1
                 continue;
             end
             temp = P(i, :);
             P(i, :) = p(i+Index(1), :);
             p(i+Index(1), :) = temp;
         end
      a=P(i,i);  
      P(i,:)= P(i,:)/a;
     for j=i+1:row     
       P(j,:)= P(j,:)- P(j,i)* P(i,:);
     end
   end
  a=P(row,row);  
  P(row,:)= P(row,:)/a;
 for i=row:-1:2   % Finding zeros of the upper triangular matrix.
     for j=i-1:-1:1    
       P(j,:)= P(j,:)- P(j,i)* P(i,:);
     end
 end 
 disp('The required solution is:')
 P
 %P(:,col)