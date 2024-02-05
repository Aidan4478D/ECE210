% tips for next homework (assignment 3)

%each dimension can be thought of a for loop
randi(bound, dim, iterations)

dice = randi(sides, rolls, iterations)
sum(dice) % can be seen as like two dice rolls


% lessons on functions!

% avoid for loops at all costs because they're slow af lol
for i = 1:10
    if i == 1
        display(i);
    end
    else
        display("poop")
    end
end

% lambda function
% anonymous function; like a function in a variable; kinda like a single line function
sq = @(x) x^2;

a = (x, y) x + y;
b = @(x) a(x, 3); 

trig = {@sin, @cos, @tan};
trig{1}(pi /6);  %will return sin(pi / 6)

% make seperate files for functions
function x = fib()

    if n <= 1
        x = 1;
        return
    end

    x = fib(n - 1) + fib(n - 2);
end



% distance formula functions
function s = sq(x)
    s = x .* x;
end

function d = dist(x, y)
    d = sqrt(sum(sq(x) + sq(y))); 
end


% outer product to yield matrix
function[sz, A] = outer_product(u v)
    A = u' * v;
    sz = size(A)
end


% look at varargin and varargout to get variable input functions
% there's documentation on classes if interested in oop stuff
