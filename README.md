# tensor2matrix
Generate matrices from tensors (Sympy indexed expressions) using Voigt and vectorization rules. 

Examples
--------
Tensor Contraction

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

	mat = generate_matrix(exprs)


Example using voigt notation, kinematic rule and vectorization
	
	B = IndexedBase("B")	# Shape function tensor (B Matrix)
	N = IndexedBase("N")	# Shape functions

	i = Idx("i", 2)			# 2D problem
    j = Idx("j", 2)
    k = Idx("k", 2)

	I = Idx("I",4)			# 4 Nodes

	expr = Equals(B[i,j,r,I], Rational(1,2)*(N[I,j]*delta(r,i) + N[I,i]*delta(r,j)))

	exprs = transform(expr, (i,j,r,I), voigt=(i,j), kinematic=(i,j), vector=(r,I))
	
	mat = generate_matrix(exprs)

The resultant matrix 

	B = [[N[0, 0],       0, N[1, 0],       0, N[2, 0],       0, N[3, 0],       0],
         [      0, N[0, 1],       0, N[1, 1],       0, N[2, 1],       0, N[3, 1]],
         [N[0, 1], N[0, 0], N[1, 1], N[1, 0], N[2, 1], N[2, 0], N[3, 1], N[3, 0]]]
