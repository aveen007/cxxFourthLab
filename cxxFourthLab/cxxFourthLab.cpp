// cxxFourthLab.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "sparseVector.h"
#include "Matrix.h"
#include <vector>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <vector>
using namespace std;

// Function to create a random dense matrix
std::vector<std::vector<int>> createDenseMatrix(int rows, int cols, float sparsity) {
    std::vector<std::vector<int>> matrix(rows, std::vector<int>(cols, 0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (static_cast<float>(rand()) / RAND_MAX > sparsity) { // Fill with non-zero values based on sparsity
                matrix[i][j] = rand() % 10 + 1; // Random values from 1 to 10
            }
        }
    }
    return matrix;
}

// Function to measure the time for an operation on a dense matrix
template<typename Func>
double measureDenseMatrixOperation(const std::vector<std::vector<int>>& matrixA,
    const std::vector<std::vector<int>>& matrixB,
    Func operation) {
    auto start = std::chrono::high_resolution_clock::now();
    operation(matrixA, matrixB);
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count(); // Return time in milliseconds
}

// Function to measure the time for an operation on a sparse matrix
template<typename Func>
double measureSparseMatrixOperation(const SparseMatrix<int>& sparseMatrixA,
    const SparseMatrix<int>& sparseMatrixB,
    Func operation) {
    auto start = std::chrono::high_resolution_clock::now();
    operation(sparseMatrixA, sparseMatrixB);
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count(); // Return time in milliseconds
}

// Example of addition operation for dense matrix
std::vector<std::vector<int>> addDenseMatrices(const std::vector<std::vector<int>>& A,
    const std::vector<std::vector<int>>& B) {
    int rows = A.size();
    int cols = A[0].size();
    std::vector<std::vector<int>> result(rows, std::vector<int>(cols, 0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

// Example of addition operation for sparse matrix
SparseMatrix<int> addSparseMatrices(const SparseMatrix<int>& A, const SparseMatrix<int>& B) {
    return A + B; // Use your overloaded + operator
}

std::vector<std::vector<int>> multiplyDenseMatrices(const std::vector<std::vector<int>>& A,
    const std::vector<std::vector<int>>& B) {
    int rowsA = A.size();
    int colsA = A[0].size();
    int colsB = B[0].size();

    // Ensure the number of columns in A equals the number of rows in B
    if (colsA != B.size()) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
    }

    std::vector<std::vector<int>> result(rowsA, std::vector<int>(colsB, 0));
    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < colsB; ++j) {
            for (int k = 0; k < colsA; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}
SparseMatrix<int> multiplySparseMatrices(const SparseMatrix<int>& A, const SparseMatrix<int>& B) {
    return A * B; 
}
std::vector<std::vector<int>> multiplyDenseMatrixByScalar(const std::vector<std::vector<int>>& A, const std::vector<std::vector<int>>& B) {
    int scalar = rand();
    int rows = A.size();
    int cols = A[0].size();

    std::vector<std::vector<int>> result(rows, std::vector<int>(cols, 0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = A[i][j] * scalar;
        }
    }
    return result;
}
SparseMatrix<int> multiplySparseMatrixByScalar(const SparseMatrix<int>& A, const SparseMatrix<int>&  B) {
    int scalar = rand();
    return A * scalar;
}

// Function to create a random sparse vector
vector<int> createDenseVector(int size, float sparsity) {
    vector<int> veector(size);
    for (int i = 0; i < size; ++i) {
            if (static_cast<float>(rand()) / RAND_MAX > sparsity) { // Fill with non-zero values based on sparsity
                veector[i] = rand() % 10 + 1; // Random values from 1 to 10
            }
        }
    return veector;
}

// Function to measure the time for an operation on a sparse vector
template<typename Func>
double measureSparseVectorOperation(const SparseVector<int>& vectorA,
    const SparseVector<int>& vectorB,
    Func operation) {
    auto start = std::chrono::high_resolution_clock::now();
    operation(vectorA, vectorB);
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count(); // Return time in milliseconds
}
template<typename Func>
double measureDenseVectorOperation(const vector<int>& vectorA,
    const vector<int>& vectorB,
    Func operation) {
    auto start = std::chrono::high_resolution_clock::now();
    operation(vectorA, vectorB);
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count(); // Return time in milliseconds
}

// Example of addition operation for sparse vector
SparseVector<int> addSparseVectors(const SparseVector<int>& A, const SparseVector<int>& B) {
    return A + B; // Use your overloaded + operator
}

// Example of dot product operation for sparse vector
int dotProductSparseVectors(const SparseVector<int>& A, const SparseVector<int>& B) {
    return A.dot(B); // Use your dot product function
}

// Example of scalar multiplication for sparse vector
SparseVector<int> multiplySparseVectorByScalar(const SparseVector<int>& A, const SparseVector<int>& B) {
    int scalar = rand();
    return A * scalar; // Use your overloaded * operator
}
std::vector<int> addVectors(const std::vector<int>& A, const std::vector<int>& B) {
    // Assuming A and B are of the same size
    std::vector<int> result(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        result[i] = A[i] + B[i];
    }
    return result;
}
int dotProductVectors(const std::vector<int>& A, const std::vector<int>& B) {
    // Assuming A and B are of the same size
    int result = 0;
    for (size_t i = 0; i < A.size(); ++i) {
        result += A[i] * B[i];
    }
    return result;
}
std::vector<int> multiplyVectorByScalar(const std::vector<int>& A, const std::vector<int>& B) {
    int scalar = rand();
    std::vector<int> result(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        result[i] = A[i] * scalar;
    }
    return result;
}









int main() {
    const int rows = 1000;
    const int cols = 1000;
    const float sparsity = 0.99; // 90% sparse

    // Create dense matrices
    auto denseMatrixA = createDenseMatrix(rows, cols, sparsity);
    auto denseMatrixB = createDenseMatrix(rows, cols, sparsity);

    // Create sparse matrices
    SparseMatrix<int> sparseMatrixA;
    SparseMatrix<int> sparseMatrixB;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (denseMatrixA[i][j] != 0) sparseMatrixA.set(i, j, denseMatrixA[i][j]);
            if (denseMatrixB[i][j] != 0) sparseMatrixB.set(i, j, denseMatrixB[i][j]);
        }
    }

    std::cout << "====== Matrix Testing ========\n";

    // Measure time for addition on dense matrix
    double denseAddTime = measureDenseMatrixOperation(denseMatrixA, denseMatrixB, addDenseMatrices);
    std::cout << "Dense Matrix Addition Time: " << std::fixed << std::setprecision(2) << denseAddTime << " ms ";

        // Measure time for addition on sparse matrix
        double sparseAddTime = measureSparseMatrixOperation(sparseMatrixA, sparseMatrixB, addSparseMatrices);
    std::cout << "Sparse Matrix Addition Time: " << std::fixed << std::setprecision(2) << sparseAddTime << " ms\n";

    double denseMultTime = measureDenseMatrixOperation(denseMatrixA, denseMatrixB, multiplyDenseMatrices);
    std::cout << "Dense Matrix Multiplication Time: " << std::fixed << std::setprecision(2) << denseMultTime << " ms ";

    // Measure time for addition on sparse matrix
    double sparseMultTime = measureSparseMatrixOperation(sparseMatrixA, sparseMatrixB, multiplySparseMatrixByScalar);
    std::cout << "Sparse Matrix Multiplication Time: " << std::fixed << std::setprecision(2) << sparseMultTime << " ms\n";

    double denseMultSTime = measureDenseMatrixOperation(denseMatrixA, denseMatrixB, multiplyDenseMatrixByScalar);
    std::cout << "Dense Matrix Multiplication by scalar Time: " << std::fixed << std::setprecision(2) << denseMultTime << " ms ";

    // Measure time for addition on sparse matrix
    double sparseMultSTime = measureSparseMatrixOperation(sparseMatrixA, sparseMatrixB, multiplySparseMatrices);
    std::cout << "Sparse Matrix Multiplication by scalar Time: " << std::fixed << std::setprecision(2) << sparseMultTime << " ms\n";
    std::cout << "====== Vector Testing ========\n";
    int size = 10000;
    vector<int> denseVA = createDenseVector(size, sparsity);
    vector<int> denseVB = createDenseVector(size, sparsity);

    SparseVector<int> sparseVectorA (size);
    SparseVector<int> sparseVectorB(size);
    for (int i = 0; i < size; ++i) {
            if (denseVA[i] != 0) sparseVectorA.set(i, denseVA[i]);
            if (denseVB[i] != 0) sparseVectorB.set(i,  denseVB[i]);
        
    }


    // Measure time for addition on sparse vector
    double DenseVAddTime = measureDenseVectorOperation(denseVA, denseVB, addVectors);
    std::cout << "Dense Vector Addition Time: " << std::fixed << std::setprecision(2) << DenseVAddTime << " ms  ";

    double sparseVAddTime = measureSparseVectorOperation(sparseVectorA, sparseVectorB, addSparseVectors);
    std::cout << "Sparse Vector Addition Time: " << std::fixed << std::setprecision(2) << sparseVAddTime << " ms \n ";

        // Measure time for dot product on sparse vector
    double DenseVMultTime = measureDenseVectorOperation(denseVA, denseVB, dotProductVectors);
    std::cout << "Dense Vector dot product Time: " << std::fixed << std::setprecision(2) << DenseVMultTime << " ms  ";

        double dotProductTime = measureSparseVectorOperation(sparseVectorA, sparseVectorB, dotProductSparseVectors);
    std::cout << "Sparse Vector Dot Product Time: " << std::fixed << std::setprecision(2) << dotProductTime << " ms\n        ";

        // Measure time for scalar multiplication on sparse vector
        int scalar = rand() % 10 + 1; // Random scalar
        double DenseVSTime = measureDenseVectorOperation(denseVA, denseVB, multiplyVectorByScalar);
        std::cout << "Dense Vector mult scalar  Time: " << std::fixed << std::setprecision(2) << DenseVSTime << " ms  ";

    double scalarMultTime = measureSparseVectorOperation(sparseVectorA, sparseVectorB, multiplySparseVectorByScalar);
    std::cout << "Sparse Vector Multiplication by Scalar Time: " << std::fixed << std::setprecision(2) << scalarMultTime << " ms \n";


        return 0;
}



// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
