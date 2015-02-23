import sympy
import collections
from IPython.display import display
from sympy import *
from sympy.tensor import IndexedBase, Idx, Indexed
from sympy.physics.secondquant import KroneckerDelta as delta
from sympy.core.relational import Relational

init_printing()


class Equals(Relational): 
    """
    Non-functional equals for display purpose only 
    """
    rel_op = '=='


def apply_recursive(f, expr):
    """
    Applys a function to each sub-tree of the expression tree, 
    returns the new expression.
    """
    if expr.args:
        new_args = [apply_recursive(f,e) for e in expr.args]
        return f(expr.func(*new_args))
    else: return f(expr)


def expand_sums(expr, indices):
    """
    Expands sums on multiplicative terms on RHS where
    indices appear exactly twice. 
    """

    if not indices: return expr
    if not isinstance(indices, collections.Iterable): indices = [indices]

    contracted = set()

    def sum_mul_term(expr, index):
        if not isinstance(expr, sympy.Mul): return expr

        count = expr.count(index)
        if (count == 2):
            contracted.add(index)
            terms = [expr.subs(index,val) for val in range(index.lower, index.upper+1)]
            expr = sympy.Add(*terms)
        return expr

    # Expand sum for every index
    rhs = expr.rhs
    for index in indices:
        rhs = apply_recursive(lambda e: sum_mul_term(e, index), rhs)

    return Equals(expr.lhs, rhs)


def enumerate_indices(expr, indices):
    """
    Generates a list of expressions from a single
    expression by substitutng values into indices
    based on their range.
    """

    # Make sure indices is a list
    if not indices: return expr
    if not isinstance(indices, collections.Iterable): indices = [indices]

    all_expr = [expr]
    indices = list(indices)        # To allow for popping below

    # For every index
    while indices:
        index = indices.pop()

        # For every previously added expression
        num_expr = len(all_expr)
        for i in range(0,num_expr):
            expr = all_expr.pop(0)
            # Substitute 
            new_exprs = [expr.subs(index, val) for val in range(index.lower, index.upper+1)]
            # Add only if doesnt exist yet
            [all_expr.append(new) for new in new_exprs if new not in all_expr]

    return all_expr


# Pair -> Parent value 
voigt_map = {1:{(0,0):0},

             2:{(0,0):0,  (0,1):2,
                (1,0):2,  (1,1):1},

             3:{(0,0):0, (0,1):5, (0,2):4,
                (1,0):5, (1,1):1, (1,2):3,
                (2,0):4, (2,1):3, (2,2):2}}    

def apply_rules(tensor, full_exprs, voigt_pairs, kinematic_pairs, vector_pairs):
    """
    Cycles through index pairs and applies voigt/vector rules if required
    """

    # Cycle through index pairs
    indices = tensor.indices
    for pos in range(len(indices)-1,-1,-1):
        pair = tuple(indices[pos:pos+2]) 

        # Determine what rule to apply
        rule = None
        rhs_transform = lambda num_pair, x: x
        if pair in voigt_pairs:
            # Get corrent dimension
            dim = pair[0].upper+1

            # Get factor to multiply rhs by
            if pair in kinematic_pairs: 
                rhs_transform = lambda num_pair, x: 2*x if num_pair[0] != num_pair[1] else x

            # Pair mapping rule
            rule = lambda num_pair: voigt_map[dim][num_pair]

        elif pair in vector_pairs:
            # Flip pair
            if (pair[0].label.name.islower()): 
                dim = pair[0].upper+1
                rule = lambda num_pair: num_pair[1]*dim + num_pair[0]
            else: 
                dim = pair[1].upper+1
                rule = lambda num_pair: num_pair[0]*dim + num_pair[1]

        # Dont transform this rule if nothing found
        if not rule: continue

        results = []
        for expr in full_exprs:
            lhs = expr.lhs
            num_pair = tuple(lhs.indices[pos:pos+2])
            new_val = rule(num_pair)
            new_indices = lhs.indices[0:pos] + (new_val,) + lhs.indices[pos+2:]
            lhs = lhs.base[new_indices]
            rhs = rhs_transform(num_pair, expr.rhs)
            results.append(Equals(lhs, rhs))
        full_exprs = results

    return full_exprs


def generate_matrix(exprs):
    """
    Converts a set of expressions to a matrix
    """

    # Count required indicies
    nindices = min([len(expr.lhs.indices) for expr in exprs])
    nrows = max([expr.lhs.indices[0] for expr in exprs]) + 1
    if nindices == 1: 
        ncols = 1
    elif nindices == 2: 
        ncols = max([expr.lhs.indices[1] for expr in exprs]) + 1
    else: 
        return exprs    # Cannot generate matrix

    #lhs_mat = zeros(nrows, ncols)
    rhs_mat = zeros(nrows, ncols)

    for expr in exprs:
        row = expr.lhs.indices[0]
        if nindices == 2: col = expr.lhs.indices[1]
        else: col = 0
        #lhs_mat[row,col] = expr.lhs
        if not rhs_mat[row,col]: rhs_mat[row,col] = expr.rhs

    return Equals(expr.lhs.base, rhs_mat)

def generate_subs_map(exprs):
    return {expr.rhs : expr.lhs for expr in exprs}

def transform(expr, indices, voigt=[], kinematic=[], vector=[]):
    """
    Converts a tensor expression to a set of expanded expressions with
    voigt and vector rules applied
    """

    # Convert voigt, kinematic, vector to lists of pairs, if required
    def ensure_list_of_pairs(x): return [x] if (x and isinstance(x[0], Idx)) else x
    voigt     = ensure_list_of_pairs(voigt)
    kinematic = ensure_list_of_pairs(kinematic)
    vector    = ensure_list_of_pairs(vector)

    # This is the lhs of the original tensor expression, which is used to determine
    # where the original indices were 
    tensor = expr.lhs
    if not isinstance(tensor, Indexed): raise TypeError("Require an indexed expression on lhs.")

    # Expand sums, generate expressions
    exprs = expand_sums(expr, indices)
    exprs = enumerate_indices(exprs, indices)
    exprs = apply_rules(tensor, exprs, voigt, kinematic, vector)
    return exprs

A = IndexedBase("A")
B = IndexedBase("B")
C = IndexedBase("\mathbb{C}")
N = IndexedBase("N")
dNdX = IndexedBase("\\frac{dN}{dX}")
IoI = IndexedBase("I \otimes I")
II = IndexedBase("\mathbb{I}")
P = IndexedBase("\mathbb{P}")
dRdu = IndexedBase("dRdu")

strain = IndexedBase("varepsilon")

i = Idx("i",2)
j = Idx("j",2)
k = Idx("k",2)
l = Idx("l",2)
r = Idx("r",2)
m = Idx("m",2)
n = Idx("n",2)

I = Idx("I",4)
J = Idx("J",4)
M = Idx("M",4)

# Strain tensor
expr = Equals(strain[i,j], strain[i,j])
exprs = transform(expr, (i,j), (i,j), (i,j))
mat = generate_matrix(exprs)
display(mat)

# B Matrix
expr = Equals(B[i,j,r,I], Rational(1,2)*(N[I,j]*delta(r,i) + N[I,i]*delta(r,j)))
exprs = transform(expr, (i,j,r,I), voigt=(i,j), kinematic=(i,j), vector=(r,I))
B = generate_matrix(exprs)
display(B)
B = B.rhs

# Trace operatpr
expr = Equals(IoI[i,j,k,l], delta(i,j)*delta(k,l))
exprs = transform(expr, (i,j,k,l), voigt=[(i,j), (k,l)])
mat = generate_matrix(exprs)
display(mat)

# 4th order identity
expr = Equals(II[i,j,k,l], 0.5*delta(i,k)*delta(j,l) + 0.5*delta(i,l)*delta(j,k))
exprs = transform(expr, (i,j,k,l), voigt=[(i,j), (k,l)])
mat = generate_matrix(exprs)
display(mat)

# Projection
expr = Equals(P[i,j,k,l], Rational(1,2)*(delta(i, k)*delta(j, l)+delta(i, l)*delta(j, k))-Rational(1,3)*delta(i, j)*delta(k, l))
exprs = transform(expr, (i,j,k,l), voigt=[(i,j), (k,l)])
mat = generate_matrix(exprs)
display(mat)

# Tangent matrix
expr = Equals(C[i,j,k,l], C[i,j,k,l])
exprs = transform(expr, (i,j,k,l), voigt=[(i,j), (k,l)])
sub_map = generate_subs_map(exprs)
CC = generate_matrix(exprs)
CC = CC.rhs

# Derivative of residual  w.r.t. displacement
expr = Equals(dRdu[I, i, J, j], N[I,k]*C[i, k, j, n]*N[J,n])
exprs = transform(expr, (i,I,j,J,k,n), vector=[(I,i), (J,j)])
exprs = [e.subs(sub_map) for e in exprs]
mat = generate_matrix(exprs)
display(mat)


