import math
from math import sqrt
import numbers

def zeroes(height, width):
        """
        Creates a matrix of zeroes.
        """
        g = [[0.0 for _ in range(width)] for __ in range(height)]
        return Matrix(g)

def identity(n):
        """
        Creates a n x n identity matrix.
        """
        I = zeroes(n, n)
        for i in range(n):
            I.g[i][i] = 1.0
        return I

class Matrix(object):

    # Constructor
    def __init__(self, grid):
        self.g = grid
        self.h = len(grid)
        self.w = len(grid[0])

    #
    # Primary matrix math methods
    #############################
        
    def determinant(self):
        """
        Calculates the determinant of a matrix.
        """
        if not self.is_square():
            raise ValueError("Cannot calculate determinant of non-square matrix.")
        
        if self.h == 1:
            return self.g[0][0]
        if self.h == 2:
            a = self.g[0][0]
            b = self.g[0][1]
            c = self.g[1][0]
            d = self.g[1][1]
            return a*d - c*b
        total = 0
        flag = False
        for l in range(self.w):
            factor = self.g[0][l]
            smallerdet = []
            for i in range(1, self.h):
                row = []
                for j in range(self.w):
                    if j != l:
                        row.append(self.g[i][j])
                smallerdet.append(row)
            det = Matrix(smallerdet)
            if flag:
                total -= factor * det.determinant()
                flag = False
            else:
                total += factor * det.determinant()
                flag = True
        return total
                            
    def deepCopy(self):
        """Returns a fresh copy of the Matrix class
        """
        newCopy = []
        for i in range(self.h):
            row = []
            for j in range(self.w):
                row.append(self.g[i][j])
            newCopy.append(row)
        return Matrix(newCopy)
    
    def __minor(self):
        """
        Creates the minor of the given 3x3 matrix and returns a fresh copy
        """
        copyMtrx = self.deepCopy()
        for l in range(self.h):
            for k in range(self.w):
                newMtrx = []
                for i in range(self.h):
                    if l != i:
                        row = []
                        for j in range(self.w):
                            if k != j:
                                row.append(self.g[i][j])
                        newMtrx.append(row)
                subMtrx = Matrix(newMtrx)
                det = subMtrx.determinant()
                copyMtrx.g[l][k] = det
        return copyMtrx
    
    def __cofactor(self):
        """
        Creates a matrix of cofactors and returns a fresh copy
        """
        cpyMtrx = self.deepCopy()
        flag = False
        for i in range(self.h):
            for j in range(self.w):
                if flag:
                    cpyMtrx.g[i][j] = - cpyMtrx.g[i][j]
                    flag = False
                else:
                    flag = True
        return cpyMtrx
    
    def trace(self):
        """
        Calculates the trace of a matrix (sum of diagonal entries).
        """
        if not self.is_square():
            raise ValueError("Cannot calculate the trace of a non-square matrix.")

        total = 0
        for i in range(self.h):
            total += self.g[i][i]
        return total
    
    def inverse(self):
        """
        Calculates the inverse of a 1x1 or 2x2 Matrix.
        """
        if not self.is_square():
            raise ValueError("Non-square Matrix does not have an inverse.")
        if self.h == 1:
            return Matrix([[1/self.g[0][0]]])
        if self.h == 2:
            det = self.determinant()
            if det == 0:
                raise ValueError("This matrix is not invertable")
            else:
                factor = 1/det
                trc = self.trace()
                identMtrx = identity(self.w)
                inverseMatrix = []
                for i in range(self.h):
                    row = []
                    for j in range(self.w):
                        row.append(factor*(trc*identMtrx[i][j] - self.g[i][j]))
                    inverseMatrix.append(row)
                return Matrix(inverseMatrix)
        # Math thanks to http://wwwf.imperial.ac.uk/metric/metric_public/matrices/inverses/inverses2.html 
        #Step 1: Compute the minor:
        inverseMtrx = self.__minor()
        #Step 2: Compute the cofactor:
        inverseMtrx = inverseMtrx.__cofactor()
        #Step 3: Transpose:
        inverseMtrx = inverseMtrx.T()
        #Step 4: Obtain determinant of the original matrix
        det = self.determinant()
        #Step 5: divide the original matrix by the determinant
        if det == 0:
            raise ValueError("This matrix is not invertable")
        else:
            inverseMtrx = 1/det * inverseMtrx
            return inverseMtrx
                
    def T(self):
        """
        Returns a transposed copy of this Matrix.
        """
        transposed = []
        for j in range(self.w):
            row = []
            for i in range(self.h):
                row.append(self.g[i][j])
            transposed.append(row)
        return Matrix(transposed)
                           
    def is_square(self):
        return self.h == self.w

    #
    # Begin Operator Overloading
    ############################
    def __getitem__(self,idx):
        """
        Defines the behavior of using square brackets [] on instances
        of this class.

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > my_matrix[0]
          [1, 2]

        > my_matrix[0][0]
          1
        """
        return self.g[idx]

    def __repr__(self):
        """
        Defines the behavior of calling print on an instance of this class.
        """
        s = ""
        for row in self.g:
            s += " ".join(["{} ".format(x) for x in row])
            s += "\n"
        return s

    def __add__(self,other):
        """
        Defines the behavior of the + operator
        """
        if self.h != other.h or self.w != other.w:
            raise ValueError("Matrices can only be added if the dimensions are the same")
            
        resultMtrx = []
        for i in range(self.h):
            row = []
            for j in range(self.w):
                row.append(self.g[i][j] + other.g[i][j])
            resultMtrx.append(row)
        return Matrix(resultMtrx)

    def __neg__(self):
        """
        Defines the behavior of - operator (NOT subtraction)

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > negative  = -my_matrix
        > print(negative)
          -1.0  -2.0
          -3.0  -4.0
        """
        resultMtrx = []
        for i in range(self.h):
            row = []
            for j in range(self.w):
                row.append(-self.g[i][j])
            resultMtrx.append(row)
        return Matrix(resultMtrx)

    def __sub__(self, other):
        """
        Defines the behavior of - operator (as subtraction)
        """
        resultMtrx = []
        for i in range(self.h):
            row = []
            for j in range(self.w):
                row.append(self.g[i][j] - other.g[i][j])
            resultMtrx.append(row)
        return Matrix(resultMtrx)

    def __mul__(self, other):
        """
        Defines the behavior of * operator (matrix multiplication)
        """
        otherTranspose = other.T()
        resultMtrx = []
        for i in range(self.h):
            row = []
            for j in range(otherTranspose.h):
                row1 = self.g[i]
                row2 = otherTranspose[j]
                total = 0
                for l in range(len(row1)):
                    total += row1[l] * row2[l]
                row.append(total)
            resultMtrx.append(row)
        return Matrix(resultMtrx)

    def __rmul__(self, other):
        """
        Called when the thing on the left of the * is not a matrix.

        Example:

        > identity = Matrix([ [1,0], [0,1] ])
        > doubled  = 2 * identity
        > print(doubled)
          2.0  0.0
          0.0  2.0
        """
        if isinstance(other, numbers.Number):
            resultMtrx = []
            for i in range(self.h):
                row = []
                for j in range(self.w):
                    row.append(other*self.g[i][j])
                resultMtrx.append(row)
            return Matrix(resultMtrx)
        else:
            raise TypeError("Expected a scaler on the left and a Matrix on the right side of the * operator.")
            