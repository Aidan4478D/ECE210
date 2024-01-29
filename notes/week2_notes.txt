clc; 
clear; 

% how to generate a numerical range
x = [0, 1, 2, 3, 4];
x = 0:4;

% two colons, step is in the middle
y = 2:-2:-12;

% want to generate vector between 0 and 1 with 20 samples
t = linspace(0, 1, 20);

% element indexing starts at 1
x(3);            % will index the third element in the array
x([1, 3, 6]);    % will return first, third, and six
x(end);          % get the end index

A = [1, 5 
     6, 2];

A(1, 2); % will return 5
A(1, :); % will get the first column

% getting properties of atrices

% gets a uniformly distributed matrix of integers (this will be 2x4 matrix)
B = randi(4, 2, 4);

length(B); % will get the largest dimension
size(B);   % will get the rows x columns
numel(B);  % get number of elements in the matrix

% what if we want a specific sized matrix tho
A = [1, 2, 3; 4, 5, 6; 7, 8, 9];
a = 1:9;

A = reshape(a, [3, 3]);  % these two functions are the same
A = reshape(a, 3, []);

a = 1:8;
A = reshape(a, 2, 2, []); % will produce 3 dimensional matrix

% elements stored in column matrix order

A = reshape(a, [3, 3])'; % can transpose the matrix and will print in linear order now


% broadcasting, specifically has to be a row column and a column vector
u = [1, 2];
v = [3; 4; 5];

u + v; % same as [u; u; u] + [v, v]

% on homework maybe use broadcasting

sum(A);        % will sum the columns down
sum(A, 'all'); % will sum all of the elements together
sum(A(1:end)); % will get the same answer but is inefficient


% will print out the time elapsed for the function to run
tic; sum(A(1:end)); toc

% will produce a 4x4 identity matrixw
eye(4);

% will produce a 4x4 zero matrix
zeros(4);
zeros(1, 2); % will create a 1x2 zero matrix

% 4x4 matrix filled with ones
ones(4);

% is a lambda function
twos = @(x) ones(x) + 1; 
