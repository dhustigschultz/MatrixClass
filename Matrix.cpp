#include "Matrix.h"


// Determines whether the matrix is square
bool Matrix::isSquare() const {
	return m_rowSize == m_colSize;
}

// Determines whether the matrix is orthogonal
bool Matrix::isOrthogonal() const {
	return transpose() == inverse();
}

// The operator() will allow easy accessing of the elements of a Matrix
double& Matrix::operator()(const std::size_t& rowNum, const std::size_t& colNum) {
	return m_elements[rowNum][colNum];
}

bool Matrix::operator==(const Matrix& other) const {
	return ((m_rowSize == other.m_rowSize) && (m_colSize == other.m_colSize) && (m_elements == other.m_elements));
}

// Arithmetic operators between a Matrix and a constant
Matrix& Matrix::operator+=(double val) {
	for (std::size_t i{ 0 }; i < m_colSize; ++i)
		for (std::size_t j{ 0 }; j < m_rowSize; ++j)
			m_elements[i][j] += val;

	return *this;
}

Matrix Matrix::operator+(double val) {
	Matrix sum{ *this }; // Not sure if this will work, but try it first. 
	// If not, then change this to use the second constructor in Matrix.h!
	sum += val;
	return sum;
}

Matrix& Matrix::operator-=(double val) {
	for (std::size_t i{ 0 }; i < m_colSize; ++i)
		for (std::size_t j{ 0 }; j < m_rowSize; ++j)
			m_elements[i][j] -= val;

	return *this;
}

Matrix Matrix::operator-(double val) {
	Matrix diff{ *this };
	diff -= val;
	return diff;
}

Matrix& Matrix::operator*=(double val) {
	for (std::size_t i{ 0 }; i < m_colSize; ++i)
		for (std::size_t j{ 0 }; j < m_rowSize; ++j)
			m_elements[i][j] *= val;

	return *this;
}

Matrix Matrix::operator*(double val) {
	Matrix product{ *this };
	product *= val;
	return product;
}

Matrix& Matrix::operator/=(double val) {
	if (val == 0)
		throw std::runtime_error("Cannot divide by zero!");

	for (std::size_t i{ 0 }; i < m_colSize; ++i)
		for (std::size_t j{ 0 }; j < m_rowSize; ++j)
			m_elements[i][j] /= val;

	return *this;
}

Matrix Matrix::operator/(double val) {
	if (val == 0)
		throw std::runtime_error("Cannot divide by zero!");

	Matrix div{ *this };
	div /= val;
	return div;
}

// Arithmetic operators between two Matrix objects
Matrix& Matrix::operator+=(const Matrix& other) {

	for (std::size_t i{ 0 }; i < m_colSize; ++i)
		for (std::size_t j{ 0 }; j < m_rowSize; ++j)
			m_elements[i][j] += other.m_elements[i][j];

	return *this;
}

Matrix Matrix::operator+(const Matrix& other) {
	Matrix sum{ *this };
	sum += other;
	return sum;
}

Matrix& Matrix::operator-=(const Matrix& other) {

	for (std::size_t i{ 0 }; i < m_colSize; ++i)
		for (std::size_t j{ 0 }; j < m_rowSize; ++j)
			m_elements[i][j] -= other.m_elements[i][j];

	return *this;
}

Matrix Matrix::operator-(const Matrix& other) {
	Matrix diff{ *this };
	diff -= other;
	return diff;
}

Matrix Matrix::operator*(const Matrix& other) {
	// Make sure the inner matrix dimensions agree
	if (m_colSize != other.m_rowSize)
		throw std::runtime_error("Inner matrix dimensions do not agree.");

	// Need to figure out what the new dimensions will be, based on 
	// the outer dimensions of the two Matrix objects being multiplied
	std::size_t prodRowSize = m_rowSize;
	std::size_t prodColSize = other.m_rowSize;

	// Then use the first constructor, to make a new blank matrix of that size
	Matrix prod(prodRowSize, prodColSize, 0.0);

	// Then do the muliplication!
	for (std::size_t i{ 0 }; i < m_rowSize; ++i) 
		for (std::size_t j{ 0 }; j < other.m_colSize; ++j) 
			for (std::size_t k{ 0 }; k < other.m_rowSize; ++k) 
				prod.m_elements[i][j] += m_elements[i][k] * other.m_elements[k][j];

	return prod;	
}


// Access functions
std::size_t Matrix::getRowSize() const {
	return m_rowSize;
}

std::size_t Matrix::getColSize() const {
	return m_colSize;
}

// Matrix-specific operations
Matrix Matrix::transpose() const {

	Matrix Transpose(m_colSize, m_rowSize, 0.0f);

	for (std::size_t i{ 0 }; i < m_colSize; ++i)
		for (std::size_t j{ 0 }; j < m_rowSize; ++j)
			Transpose(i, j) = m_elements[j][i];

	return Transpose;
}

double Matrix::determinant() const {
	// Check to make sure that we have a square matrix first
	if (!isSquare())
		throw std::runtime_error("Matrix is not square.");

	// If there are no rows/columns in Matrix, return -1
	if (m_rowSize == 0)
		return -1.0; // Note: try to use std::optional instead, to return something indicating that it doesn't have a determinant.

	// If we have only one element in Matrix, return that element.
	if (m_rowSize == 1)
		return m_elements[0][0];

	// Formula for a 2x2 Matrix
	if (m_rowSize == 2)
		return m_elements[0][0] * m_elements[1][1] - m_elements[0][1] * m_elements[1][0];

	// Matrix is larger than 2x2:
	double det{ 0 };
	int sign{ 1 };
	for (std::size_t i{ 0 }; i < m_rowSize; ++i) {

		// First, take a submatrix of m_elements
		Matrix subMatrix(m_rowSize - 1, m_colSize - 1, 0.0f);
		for (std::size_t m{ 1 }; m < m_rowSize; ++m) {
			std::size_t z{ 0 };
			for (std::size_t n{ 0 }; n < m_colSize; ++n) {
				if (n != i) {
					subMatrix(m - 1, z) = m_elements[m][n];
					++z;
				}
			}
		}

		// Then, make a recursive call, 
		// and use the results to calculate the determinant
		det = det + sign * m_elements[0][i] * subMatrix.determinant();
		sign = -sign;
	}

	return det;
}

Matrix Matrix::cofactor() const {
	if (!isSquare())
		throw std::runtime_error("Matrix is not square.");

	Matrix solution(m_rowSize, m_colSize, 0.0f);
	Matrix subMatrix(m_rowSize - 1, m_colSize - 1, 0.0f);

	for (std::size_t i{ 0 }; i < m_rowSize; ++i) {

		for (std::size_t j{ 0 }; j < m_colSize; ++j) {

			std::size_t p{ 0 };
			for (std::size_t x{ 0 }; x < m_rowSize; ++x) {
				if (x == i)
					continue;
				std::size_t q{ 0 };

				for (std::size_t y{ 0 }; y < m_rowSize; ++y) {
					if (y == j)
						continue;

					subMatrix(p, q) = m_elements[x][y];
					++q;
				}
				++p;
			}
			solution(i, j) = static_cast<double>(pow(-1, i + j)) * subMatrix.determinant();
		}
	}
	return solution;
}

Matrix Matrix::inverse() const {
	if (determinant() == 0.0)
		throw std::runtime_error("Determinant is zero.");

	double d{ 1.0 / (determinant()) };
	Matrix matCopy(m_rowSize, m_colSize, 0.0);
	Matrix solution(m_rowSize, m_colSize, m_elements);

	solution = solution.cofactor();
	solution = solution.transpose();

	for (std::size_t i{ 0 }; i < m_rowSize; ++i)
		for (std::size_t j{ 0 }; j < m_colSize; ++j) {
			solution(i, j) = solution(i, j) * d;

			// I seem to be getting negative zeros in my resulting inverse matrix.
			// Signed zero is supported in IEEE floating point standard, and somehow 
			// that's happening somewhere in the above functions. Not sure yet why. 
			// So, for now, I'm just going to take care of that issue here.
			if (solution(i, j) == -0.0)
				solution(i, j) = -1 * solution(i, j);
		}

	return solution;
}

void Matrix::print() const {

	for (std::size_t i{ 0 }; i < m_elements.size(); ++i) {
		for (std::size_t j{ 0 }; j < m_elements[i].size(); ++j)
			std::cout << m_elements[i][j] << ' ';
		std::cout << '\n';
	}
	std::cout << '\n' << '\n';
}