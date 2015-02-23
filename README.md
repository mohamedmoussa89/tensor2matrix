# tensor2mat
Generate matrices from tensors (Sympy indexed expressions) using Voigt and vectorization rules. 

Example
-------
    A = IndexedBase("A")
    B = IndexedBase("B")
    C = IndexedBase("C")
    
    i = Idx("i", 2)
    j = Idx("j", 2)
    k = Idx("k", 2)

    expr = Equals(C[i,j], A[i,k]*B[k,j])
    
    exprs = transform(expr, (i,j,k))

The result is then

    C[0, 0] == A[0, 0]*B[0, 0] + A[0, 1]*B[1, 0]
	C[1, 0] == A[1, 0]*B[0, 0] + A[1, 1]*B[1, 0]
	C[0, 1] == A[0, 0]*B[0, 1] + A[0, 1]*B[1, 1]
	C[1, 1] == A[1, 0]*B[0, 1] + A[1, 1]*B[1, 1]

This can be converted to a matrix using

	generate_matrix(exprs)
