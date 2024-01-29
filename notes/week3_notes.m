x = 1:8;
A = reshape(1:16, [4,4]);


x([1, 1, 1]); % will return 1 1 1
x(3:5) = 0;   % sets the 3, 4, 5 index = 0


A(3, 3); % get the elemen in third row and third column


A([1,4], [2,3]); % get elements in rows 1 and 4 and columns 2 and 3
A([1,4], 1:4);   % get elements in rows 1 and 4 and all columns


x(1:end); % get a row vector of all elements
x(:);     % get a column vector of all elements


A(:, 1);  % get all elements from column 1


% has column major ordering; indexes down columns not across rows
A(11); % will return 11

% get the index (row, col) of a specific element
[row, col] = ind2sub(size(A), 6);

% get element at a specific index
ind = sub2ind(size(A), row, col);

% any value that is non-zero is considered to be true
corners = logical([1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]);

A(corners); % will return 1, 4, 13, 16; corners in the array A

all(corners); % will return false (there are some 0 elements) 
any(corners); % will return true (there are some 1 elements)

% short circuiting = stop if anything is false

A == 9;              % returns a logical array with index 9 as 1 and everything else as 0
A == [3, 6, 11, 15];
A < 7;               % same sense; returns a logical array with every index < 7 as a 1 and everything else 0
A < 7 & A > 3;       % same thing just tighter bounds

% to get the values do like
A(A < 7 & A > 3) % will get back elements that satisfy the conditions

% repeat matrix
repmat(A, 4, 1);
