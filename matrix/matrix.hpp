/*
    Implementation of Arduino Matrix class
    @author: Joe Mattekatt
*/

#include "matrix.h"


/*
    Default constructor. Memory for the Matrix is not allocated.
*/
template <class T>
Matrix<T>::Matrix() {
    rows = DEFAULT_VALUE;
    cols = DEFAULT_VALUE;
    size = rows * cols;
    head = nullptr;
}

/*
    Constructor. Memory for the Matrix is allocated; values initialized to zero.
    @param init_rows: Number of rows in Matrix
    @param init_cols: Number of columns in Matrix
*/
template <class T>
Matrix<T>::Matrix(const unsigned int init_rows, const unsigned int init_cols) {
    rows = init_rows;
    cols = init_cols;
    size = rows * cols;

    head = new T[size];
    for (int i = START_INDEX; i < size; ++i) {
        this->at(i) = DEFAULT_VALUE;
    }
}

/*
    Constructor. Memory for the Matrix is allocated; values initialized to fill_value.
    @param init_rows: Number of rows in Matrix
    @param init_cols: Number of columns in Matrix
    @param fill_value: Value used to initialize Matrix elements
*/
template <class T>
Matrix<T>::Matrix(const unsigned int init_rows, const unsigned int init_cols,
    const T fill_value) : Matrix(init_rows, init_cols) {
    for (int i = START_INDEX; i < size; ++i) {
        this->at(i) = fill_value;
    }
}

/*
    Constructor. Memory for the Matrix is allocated; values initialized to init_vals.
    @param init_rows: Number of rows in Matrix
    @param init_cols: Number of columns in Matrix
    @param init_vals: Values used to initialize Matrix elements in row major order
*/
template <class T>
Matrix<T>::Matrix(const unsigned int init_rows, const unsigned int init_cols,
    const T* init_vals) : Matrix(init_rows, init_cols) {
    ASSERT((init_vals != nullptr), "ERROR: Initial values array cannot be nullptr!");
    for (int i = START_INDEX; i < size; ++i) {
        this->at(i) = init_vals[i];
    }
}


/*
    Destructor.
*/
template <class T>
Matrix<T>::~Matrix() {
    if (head != nullptr) {
        delete[] head;
    }
}

/*
    Helper function to copy another Matrix's members into the current Matrix.
    Leaves the other Matrix with default member values and the current Matrix's
    storage array.
    @param other: The Matrix to copy members from
*/
template <class T>
inline void Matrix<T>::swap(Matrix<T> &other) {
    rows = other.rows;
    other.rows = DEFAULT_VALUE;
    cols = other.cols;
    other.cols = DEFAULT_VALUE;
    size = other.size;
    other.size = DEFAULT_VALUE;

    T* tmp = head;
    head = other.head;
    other.head = tmp;
}

/*
    Move constructor. Moves the other Matrix's members into the current matrix.
    Leaves the other Matrix in default state.
    @param other: The Matrix to move members from
*/
template <class T>
Matrix<T>::Matrix(Matrix<T> &&other) 
    : Matrix() {
    this->swap(other);
}

/*
    Copy constructor. Copies the other Matrix's members into the current matrix.
    @param other: The Matrix to copy members from
*/
template <class T>
Matrix<T>::Matrix(const Matrix<T> &other) 
    : Matrix(other.rows, other.cols, other.head) {
}

/*
    Assignment operator. Moves the other Matrix's members into the current matrix 
    given rvalue reference to other Matrix. Copies the other Matrix's members into
    the current matrix given lvalue reference to other Matrix.
    @param other: The Matrix to move/copy members from; can be passed as rvalue/lvalue
    @return: The current Matrix
*/
template <class T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> other) {
    this->swap(other);

    return *this;
}

/*
    Print the elements of the Matrix distributed in rows and columns.
*/
template <class T>
void Matrix<T>::print() const {
    for (int i = START_INDEX; i < rows; ++i) {
        for (int j = START_INDEX; j < cols; ++j) {
            Serial.print(this->at(i, j));
            Serial.print("\t");
        }
        Serial.println();
    }
}

/*
    Element access helper function. Returns read-only element at [row, column] position.
    @param row: Element's row in Matrix; in range of (0, rows-1)
    @param col: Element's column in Matrix; in range of (0, cols-1)
    @return: Element of interest
*/
template <class T>
inline T Matrix<T>::at(const unsigned int row, const unsigned int col) const {
    ASSERT((head != nullptr), "ERROR: Matrix has no memory allocated!");
    return head[row * cols + col];
}

/*
    Element access helper function. Returns writable element at [row, column] position.
    @param row: Element's row in Matrix; in range of (0, rows-1)
    @param col: Element's column in Matrix; in range of (0, cols-1)
    @return: Reference to element of interest
*/
template <class T>
inline T& Matrix<T>::at(const unsigned int row, const unsigned int col) {
    ASSERT((head != nullptr), "ERROR: Matrix has no memory allocated!");
    return head[row * cols + col];
}

/*
    Element access helper function. Returns read-only element at index in row major form.
    @param index: Position of element in row major array; in range of (0, size-1)
    @return Element of interest
*/
template <class T>
inline T Matrix<T>::at(const unsigned int index) const {
    return head[index];
}

/*
    Element access helper function. Returns writable element at index in row major form.
    @param index: Position of element in row major array; in range of (0, size-1)
    @return Reference to element of interest
*/
template <class T>
inline T& Matrix<T>::at(const unsigned int index) {
    return head[index];
}

/*
    Helper function to apply scalar operations to Matrix.
    @param scalar: Number to add, subtract or multiply with elements in Matrix
    @param funct: Function that denotes the operation to apply
    @return: New Matrix with scalar operation applied
*/
template <class T>
template <typename U>
inline decltype(auto) Matrix<T>::apply_scalar_operator(const U scalar, auto funct) const {
    Matrix<decltype(this->at(START_INDEX)*scalar)> rtn(rows, cols);

    for (int i = START_INDEX; i < size; ++i) {
        rtn.at(i) = funct(this->at(i), scalar);
    }
    return rtn;
}

/*
    Scalar addition.
    @param scalar: Number to add to all elements in Matrix
    @return: New Matrix with scalar addition applied
*/
template <class T>
template <typename U>
decltype(auto) Matrix<T>::operator+(const U scalar) const {
    return apply_scalar_operator(scalar, [](T x, U y){return x+y;});
}

/*
    Scalar subtraction.
    @param scalar: Number to subtract from all elements in Matrix
    @return: New Matrix with scalar subtraction applied
*/
template <class T>
template <typename U>
decltype(auto) Matrix<T>::operator-(const U scalar) const {
    return apply_scalar_operator(scalar, [](T x, U y){return x-y;});
}

/*
    Scalar multiplication.
    @param scalar: Number to multiply with all elements in Matrix
    @return: New Matrix with scalar multiplication applied
*/
template <class T>
template <typename U>
decltype(auto) Matrix<T>::operator*(const U scalar) const {
    return apply_scalar_operator(scalar, [](T x, U y){return x*y;});
}

/*
    Helper function to apply matrix operations to Matrix.
    @param other: Matrix to add, subtract or multiply with current Matrix
    @param funct: Function that denotes the operation to apply
    @return: New Matrix with matrix operation applied
*/
template <class T>
template <typename U>
inline decltype(auto) Matrix<T>::apply_matrix_operator(const Matrix<U> &other,
    auto funct) const {
    Matrix<decltype(this->at(START_INDEX)*other.at(START_INDEX))> rtn(rows, cols);

    for (int i = START_INDEX; i < size; ++i) {
        rtn.at(i) = funct(this->at(i), other.at(i));
    }
    return rtn;
}

/*
    Matrix addition.
    @param other: Matrix to add to current Matrix
    @return: New Matrix with matrix addition applied
*/
template <class T>
template <typename U>
decltype(auto) Matrix<T>::operator+(const Matrix<U> &other) const {
    ASSERT((rows == other.rows && cols == other.cols),
        "ERROR: Matrices must be the same size to add!");
    return apply_matrix_operator(other, [](T x, T y){return x+y;});
}

/*
    Matrix subtraction.
    @param other: Matrix to subtract from current Matrix
    @return: New Matrix with matrix subtraction applied
*/
template <class T>
template <typename U>
decltype(auto) Matrix<T>::operator-(const Matrix<U> &other) const {
    ASSERT((rows == other.rows && cols == other.cols), 
        "ERROR: Matrices must be the same size to subtract!");
    return apply_matrix_operator(other, [](T x, T y){return x-y;});
}

/*
    Matrix multiplication.
    @param other: Matrix to multiply with current Matrix
    @return: New Matrix with matrix multiplication applied
*/
template <class T>
template <typename U>
decltype(auto) Matrix<T>::operator*(const Matrix<U> &other) const {
    ASSERT((cols == other.rows), ("ERROR: Number of columns in first matrix must"
        " equal number of rows in second matrix to multiply!"));

    Matrix<decltype(this->at(START_INDEX)*other.at(START_INDEX))> rtn(rows, other.cols);

    for (int i = START_INDEX; i < rows; ++i) {
        for (int j = START_INDEX; j < other.cols; ++j) {
            
            for (int k = START_INDEX; k < cols; ++k) {
                rtn.at(i, j) += this->at(i, k) * other.at(k, j);
            }

        }
    }

    return rtn;
}


/*
    Matrix equality comparison.
    @param other: Matrix to compare to the current Matrix
    @return: True only if both matrices are the same shape and have elements in 
    the same positions; false otherwise.
*/
template <class T>
bool Matrix<T>::operator==(const Matrix<T> &other) const {
    ASSERT((rows == other.rows && cols == other.cols), 
        "ERROR: Matrices must be the same size to compare!");

    for (int i = START_INDEX; i < size; ++i) {
        if (this->at(i) != other.at(i)) {
            return false;
        }
    }
    
    return true;
}

/*
    Matrix transpose. Returns new Matrix with row and columns of elements swapped.
    @return: New transposed Matrix
*/
template <class T>
Matrix<T> Matrix<T>::transpose() const{
    Matrix<T> rtn(cols, rows);

    for (int i = START_INDEX; i < cols; ++i) {
        for (int j = START_INDEX; j < rows; ++j) {
            rtn.at(i, j) = this->at(j, i);
        }
    }

    return rtn;
}

/*
    Matrix LU decomposition.
    Given Matrix A, finds Matrices L and U such that L*U = A.
    L is a lower-triangular Matrix with ones in the diagonal.
    U is an upper-triangular Matrix.
    @param L: Reference to Matrix in which to store the lower-triangular Matrix
    @param U: Reference to Matrix in which to store the upper-triangular Matrix
*/
template <class T>
void Matrix<T>::LU_decompose(Matrix<float> &L, Matrix<float> &U) const {
    ASSERT((cols == rows), "ERROR: Matrix must be square to LU decompose!");

    L = static_cast<Matrix<float> &&>(Matrix<float>(rows, cols));
    U = static_cast<Matrix<float> &&>(Matrix<float>(rows, cols));

    for (int i = START_INDEX; i < rows; ++i) {
        L.at(i, i) = 1.0;
    }
    for (int j = START_INDEX; j < cols; ++j) {
        U.at(START_INDEX, j) = this->at(START_INDEX, j);
    }

    for (int i = 1; i < rows; ++i) {
        for (int j = START_INDEX; j < i; ++j) {
            float sum = DEFAULT_VALUE;
            for (int k = START_INDEX; k < j; ++k) {
                sum += U.at(k, j)*L.at(i, k);
            }
            L.at(i, j) = (this->at(i, j) - sum)/U.at(j, j);
        }

        for (int j = i; j < cols; ++j) {
            float sum = DEFAULT_VALUE;
            for (int k = START_INDEX; k < i; ++k) {
                sum += U.at(k, j)*L.at(i, k);
            }
            U.at(i, j) = this->at(i, j) - sum;
        }
    }
}

/*
    Solves Matrix system for the given solution set.
    Given Matrix A, column vector x and solution set B, finds x such that A*x = B.
    @param B: Matrix (column vector) denoting the solution set
    @return: Matrix that denotes the answer for the Matrix system
*/
template <class T>
Matrix<float> Matrix<T>::solve_for(Matrix<T> &B) const {
    ASSERT((cols == rows), "ERROR: Matrix must be square to solve!");

    Matrix<float> U(rows, cols);
    Matrix<float> L(rows, cols);
    this->LU_decompose(L, U);

    Matrix<float> y(rows, 1);
    y.at(START_INDEX) = B.at(START_INDEX);

    for (int i = 1; i < rows; ++i) {
        y.at(i) = B.at(i);
        for (int j = START_INDEX; j < i; ++j) {
            y.at(i) -= L.at(i, j) * y.at(j);
        }
        y.at(i) /= L.at(i, i);
    }

    Matrix<float> x(rows, 1);

    for (int i = rows-1; i >= START_INDEX; --i) {
        x.at(i) = y.at(i);

        for (int j = cols-1; j > i; --j) {
            x.at(i) -= U.at(i, j) * x.at(j);
        }
        x.at(i) /= U.at(i, i);
    }
    return x;
}
