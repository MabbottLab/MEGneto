function x = matrix_triangular_to_vector (B)

x = B(triu(true(size(B)),1));

return