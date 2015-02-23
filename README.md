# tensor2mat
Generate matrices from tensors (Sympy indexed expressions) using Voigt and vectorization rules. 

    A = IndexedBase("A")
    B = IndexedBase("B")
    C = IndexedBase("C")
    
    i = Idx("i", 2)
    j = Idx("j", 2)
    k = Idx("k", 2)

    expr = Equals(C[i,j], A[i,k]*B[k,j])
    
    exprs = transform(expr, (i,j,k))
    mat = generate_matrix(exprs)

