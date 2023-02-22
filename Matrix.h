#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <stdexcept>

class Matrix {

	// This class is limited to 2D matrices, 
	std::size_t m_rowSize{};
	std::size_t m_colSize{};
	std::vector<std::vector<double>> m_elements{};

public:
	// This constructor will initialize the Matrix with a vector of vectors of type float. 
	// 
	Matrix(const std::size_t rowSize, const std::size_t colSize, const std::vector<std::vector<double>>& elements)
		: m_rowSize(rowSize), m_colSize(colSize) {
		// Everything below ensures that m_rowSize and m_colSize will 
		// match with the resulting matrix, populated using 'elements'

		// Resize the rows and columns of m_elements, initialize with zeros
		m_elements.resize(m_rowSize);
		for (std::size_t i{ 0 }; i < m_rowSize; ++i)
			m_elements[i].resize(m_colSize, 0.0);

		// Then, populate m_elements with the contents of 'elements'. 
		// If the row and colum sizes of 'elements' are larger than 
		// those of m_elements, then the extra values get discarded.
		// If the row and colum sizes of 'elements' are smaller than
		// those of m_elements, then the remainder of m_elements is 
		// padded with zeros from the initialization above.
		std::size_t rows = (m_rowSize < elements.size()) ? m_rowSize : elements.size();
		for (std::size_t i{ 0 }; i < rows; ++i) {
			std::size_t cols = (m_colSize < elements[i].size()) ? m_colSize : elements[i].size();
			for (std::size_t j{ 0 }; j < cols; ++j)
				m_elements[i][j] = elements[i][j];
		}
	}

	// This constructor will initialize the Matrix with elements of value 'init'
	Matrix(const std::size_t rowSize, const std::size_t colSize, double init) : m_rowSize(rowSize), m_colSize(colSize) {
		m_elements.resize(m_rowSize);
		for (std::size_t i{ 0 }; i < m_elements.size(); ++i)
			m_elements[i].resize(m_colSize, init);
	}

	// Overloaded operators
	// The operator() will allow easy accessing of the elements of a Matrix
	double& operator()(const std::size_t& rowNum, const std::size_t& colNum);
	bool operator==(const Matrix& other) const;

	//Arithmetic operations between a Matrix and a constant
	Matrix& operator+=(double val);
	Matrix operator+(double val);

	Matrix& operator-=(double val);
	Matrix operator-(double val);

	Matrix& operator*=(double val);
	Matrix operator*(double val);

	Matrix& operator/=(double val);
	Matrix operator/(double val);


	//Arithmetic operations between two matrices
	Matrix& operator+=(const Matrix& other);
	Matrix operator+(const Matrix& other);

	Matrix& operator-=(const Matrix& other);
	Matrix operator-(const Matrix& other);

	Matrix operator*(const Matrix& other);

	// Determines whether the matrix is square
	bool isSquare() const;
	bool isOrthogonal() const;

	// Access functions
	std::size_t getRowSize() const;
	std::size_t getColSize() const;

	// Matrix-specific operations
	Matrix transpose() const;
	double determinant() const;
	Matrix cofactor() const;
	Matrix inverse() const;

	// Maybe add a trace function...

	// ToDo: If feeling ambitious, write a function to find eigenvalues, 
	// for matrices up to 5x5 using the Jacobi eigenvalue algorithm
	// which is here: https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm#Algorithm

	// Prints the matrix. 
	void print() const;

};

#endif //MATRIX_H_
