#pragma once
#include <iostream>
#include <unordered_map>
#include <utility>

struct pair_hash
{
    int operator()(const std::pair<int, int>& values) const noexcept
    {
        return (std::hash<int>()(values.first)
            ^ std::hash<int>()(values.second));
    }
};
template <typename T>
class SparseMatrix {
private:
    std::unordered_map<std::pair<int, int>, T, pair_hash> data;
    int rows;
    int cols;

public:
    SparseMatrix(int rows, int cols) : rows(rows), cols(cols) {}
    SparseMatrix() = default;

    // Set value at (row, col)
    void set(int row, int col, T value) {
        if (value != 0) {
            data[{row, col}] = value;  // Store non-zero values
        }
        else {
            data.erase({ row, col });//the element that was there is removed if it existed
        }
    
    }

    // Get value at (row, col)
    T get(int row, int col) const {
        auto it = data.find({ row, col });
        return (it != data.end()) ? it->second : T(0);  // Return zero if not found
    }

    // Get number of rows
    int getRows() const {
        return rows;
    }

    // Get number of columns
    int getCols() const {
        return cols;
    }


    T& operator()(int row, int col) {
        return data[{row, col}]; // Allow direct access and modification
    }

    const T& operator()(int row, int col) const {
        auto it = data.find({ row, col });
        if (it != data.end()) {
            return it->second;
        }
        return T(); // Return default value if not found
    }
    bool operator==(const SparseMatrix<T>& other) const {
        if (this->rows != other.rows || this->rows != other.cols) {
            return false;
        }

        for (const auto& pair : this->data) {
            if (other.data.find(pair.first) == other.data.end() || other(pair.first.first, pair.first.second) != pair.second) {
                return false;
            }
        }

        return true;
    }
    SparseMatrix<T> operator-(const SparseMatrix<T>& other) const {
        if (this->rows != other.rows || this->cols != other.cols) {
            throw std::invalid_argument("Matrices must have the same dimensions for subtraction.");
        }

        SparseMatrix<T> result(rows, cols);
        for (const auto& pair : data) {
            result.set(pair.first.first, pair.first.second, pair.second - other(pair.first.first, pair.first.second));
        }

        for (const auto& pair : other.data) {
            if (result.data.find(pair.first) == result.data.end()) {
                result.set(pair.first.first, pair.first.second, -pair.second);
            }
        }

        return result;
    }

    SparseMatrix<T> operator+(const SparseMatrix<T>& other) const {
        SparseMatrix<T> result;
        for (const auto& pair : data) {
            result.set(pair.first.first, pair.first.second,  pair.second+other(pair.first.first, pair.first.second));
        }
 
        return result;
    }

    SparseMatrix<T> transpose() const {
        SparseMatrix<T> result(cols, rows);
        for (const auto& pair : data) {
            result(pair.first.second, pair.first.first) = pair.second;
        }
        return result;
    }
    bool isSquare() const {
        return rows == cols;
    }

    SparseMatrix<T> inverse() const {

        int aug_rows = rows;
        int aug_cols = cols * 2;
        SparseMatrix<T> aug_matrix(aug_rows, aug_cols);
        if (!isSquare()) {
            throw std::invalid_argument("Matrix is singular!!");
        }
        for (auto element : data) {
            aug_matrix.set(element.first.first, element.first.second, element.second) ;
        }
        for (int i = 0; i < aug_rows; i++) {
          

            aug_matrix.set(i, rows + i, 1);
        }
        // Perform Gauss-Jordan elimination
        for (int i = 0; i < aug_rows; i++) {
            // Find pivot
            T pivot = aug_matrix.get(i, i);
            if (pivot == 0) {
                throw std::invalid_argument("Matrix is singular and cannot be inverted!");
            }

            // Normalize the pivot row
            for (int j = 0; j < aug_cols; j++) {
                aug_matrix.set(i, j, aug_matrix.get(i, j) / pivot);
            }

            // Eliminate the current column in all other rows
            for (int k = 0; k < aug_rows; k++) {
                if (k != i) {
                    T factor = aug_matrix.get(k, i);
                    for (int j = 0; j < aug_cols; j++) {
                        aug_matrix.set(k, j, aug_matrix.get(k, j) - factor * aug_matrix.get(i, j));
                    }
                }
            }
        }

        // Extract the right half (the inverted matrix)
        SparseMatrix<T> inverse_matrix(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = cols; j < aug_cols; j++) {
                T value = aug_matrix.get(i, j);
                if (value != 0) { // Only store non-zero values
                    inverse_matrix.set(i, j - cols, value);
                }
            }
        }

        return inverse_matrix;
    }


    


    SparseMatrix<T> operator^(const T& exponant) const {
     
        SparseMatrix<T> result(rows, cols);
        for (const auto& pair : data) {
            result.set(pair.first.first,pair.first.second, std::pow(pair.second , exponant))  ;
        }
        return result;
    }
    SparseVector<T> operator*(const SparseVector<T>& vec) const {
        if (cols != vec.getSize()) {
            throw std::invalid_argument("Matrix column count must match vector size");
        }

        SparseVector<T> result(rows);
        for (const auto& pair : data) {
            result[pair.first.first] += pair.second * vec[pair.first.second];
        }
        return result;
    }

    SparseMatrix<T> operator*(const SparseMatrix<T>& other) const {
        if (!isSquare()) {
            throw std::invalid_argument("Matrix column count of the first matrix must match the row count of the second matrix.");
        }

        SparseMatrix<T> result(rows, other.cols);

        // Iterate through non-zero elements of the first matrix
        for (const auto& pair1 : data) {
            int row1 = pair1.first.first;      // Row index of the first matrix
            int col1 = pair1.first.second;     // Column index of the first matrix
            T value1 = pair1.second;           // Value at (row1, col1)

            // Now iterate through non-zero elements of the second matrix
            for (const auto& pair2 : other.data) {
                int row2 = pair2.first.first;  // Row index of the second matrix
                int col2 = pair2.first.second; // Column index of the second matrix
                T value2 = pair2.second;       // Value at (row2, col2)

                // Check if the column of the first matrix matches the row of the second matrix
                if (col1 == row2) {
                    result(row1,col2) += value1 * value2;
                }
            }
        }
        return result;
    }

    SparseMatrix<T> operator*(const T& scalar) const {
        SparseMatrix<T> result(rows, cols);
        for (const auto& pair : data) {
            result.set(pair.first.first, pair.first.second, pair.second * scalar);
        }
        return result;
    }
    SparseMatrix<T> operator/(const T& denominator) const {
        if (denominator == T(0)) {
            throw std::invalid_argument("Cannot divide by zero.");
        }

        SparseMatrix<T> result(rows, cols);
        for (const auto& pair : data) {
            result.set(pair.first.first, pair.first.second, pair.second / denominator);
        }
        return result;
    }

    // Matrix exponentiation for integer exponents
    SparseMatrix<T> pow(int exponent) const {
        if (rows != cols) {
            throw std::invalid_argument("Matrix is not square!!");
        }
        SparseMatrix<T> result(rows, cols);
        SparseMatrix<T> base = *this;

        // Initialize result as the identity matrix
        for (int i = 0; i < rows; ++i) {
            result.set(i, i, T(1));
        }

        while (exponent) {
            if (exponent % 2 == 1) {
                result = result * base; // Implement a multiplication operator
            }
            base = base * base; // Implement a multiplication operator
            exponent /= 2;
        }
        return result;
    }
    SparseMatrix<T> pow(double exponent) const {
        if (rows != cols) {
            throw std::invalid_argument("Matrix is not square!!");
        }
        
        SparseMatrix<T> logMatrix =log(); // Compute the logarithm of the matrix
        SparseMatrix<T> result = logMatrix * exponent; // Scale the log matrix by the exponent
        return result.exp(); // Return the exponential of the result
    }
    // Matrix exponentiation for real exponents (using series expansion)
    SparseMatrix<T> exp() const {
        if (rows != cols) {
            throw std::invalid_argument("Matrix must be square for exponentiation");
        }

        SparseMatrix<T> result(rows, cols);
        SparseMatrix<T> term=(*this);
        for (int i = 0; i < rows; ++i) {
            result.set(i, i, T(1));
        }
        result = result + term;
        int n = 2;

        // Calculate e^(A) = I + A + (A^2)/2! + (A^3)/3! + ...
        while (n<100){ // Limit the series to a reasonable number of terms
          
                long long int f = factorial(n);
                term = term * (*this);
                term = term / static_cast<T>(n); 
                result = result + term;
            ++n;
        }
        return result;
    }
    double frobeniusNorm(const SparseMatrix<T>& matrix)const {
        double norm = 0.0;
        for (int i = 0; i <= matrix.rows; ++i) {
            for (int j = 0; j <= matrix.cols; ++j) {
                norm += std::pow(std::abs(matrix(i, j)), 2);
            }
        }
        return std::sqrt(norm);
    }
    // Factorial utility function
    static int factorial(int n) {
        return (n <= 1) ? 1 : n * factorial(n - 1);
    }
     SparseMatrix<T> log() const{
      if (!isSquare()) {
        throw std::invalid_argument("Matrix must be square to compute the logarithm.");
      }

      SparseMatrix<T> result(rows, cols);
      SparseMatrix<T> I(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; j++) {
            result.set(i, i, T(0));

            }
        }
      for (int i = 0; i < rows; ++i) {
          I.set(i, i, T(1));
      }
      if (*this == I) {

        return result;
      }

   
     
      if (frobeniusNorm(*this) >= 1)
        throw std::invalid_argument("norm must be less than 1.");

      SparseMatrix<T> A_minus_I = *this - I;
      SparseMatrix<T> term = A_minus_I;

      for (int n = 1; n <= 100; ++n) {
        if (n > 1) {
          term = term * A_minus_I;
        }
        if (n % 2 == 0) {
          result = result - term / static_cast<T>(n);
        } else {
          result = result + term / static_cast<T>(n);
        }
      }

      return result;
    }


    // Iterator class for SparseMatrix
    class Iterator {
    private:
        typename std::unordered_map<std::pair<int, int>, T>::iterator it;

    public:
        Iterator(typename std::unordered_map<std::pair<int, int>, T>::iterator iterator) : it(iterator) {}

        std::pair<std::pair<int, int>, T> operator*() {
            return *it;
        }

        Iterator& operator++() {
            ++it;
            return *this;
        }

        bool operator!=(const Iterator& other) const {
            return it != other.it;
        }
    };

    // Begin and end methods for iteration
    Iterator begin() {
        return Iterator(data.begin());
    }

    Iterator end() {
        return Iterator(data.end());
    }
};
