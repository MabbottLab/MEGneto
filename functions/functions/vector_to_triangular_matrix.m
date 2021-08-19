function b = vector_to_triangular_matrix (x)

% rts = roots([1 -1 -2*length(x)]);
% n = rts(rts>0);

n = (1+sqrt(1+8*length(x)))/2;

% a = triu(magic(n),1);
% a = nonzeros(a);

b = triu(ones(n),1);
b(b==1) = x;

b = b + b';
b = b + eye(n);

return