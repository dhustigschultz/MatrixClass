// MatrixPointMain.cpp : 

#include <iostream>
#include "Matrix.h"

namespace config {
    constexpr double pi{ 3.14159265358979323846 };
    constexpr std::size_t rowSize = 2;
    constexpr std::size_t colSize = 2;
}

int main()
{
    
    std::vector<std::vector<double>> E1 = {
        {1.0,0.0},{0.0,1.0}
    };
    std::vector<std::vector<double>> E2 = {
        {std::cos(config::pi / 2.0),-1 * std::sin(config::pi / 2.0)},{std::sin(config::pi / 2.0),std::cos(config::pi / 2.0)}
    };

    Matrix M1(config::rowSize, config::colSize, E1);
    Matrix M2(config::rowSize, config::colSize, 0.0);
    Matrix M3(config::rowSize, config::colSize, E2);


    std::cout << "Is the matrix square? " << std::boolalpha << M1.isSquare() << '\n';
    std::cout << "The matrix dimensions are: " << M1.getRowSize() << "x" << M1.getColSize() << '\n' << '\n';

    M1.print();

    std::cout << "The determinant of M1 is: " << M1.determinant() << '\n' << '\n';

    // Testing inversion of a Matrix, which also tests transpose, determinant, and cofactor functions
    M2 = M1.inverse();
    M2.print();

    std::cout << "Does M1 equal M2? " << (M1 == M2) << '\n' << '\n';
    std::cout << "Is M1 orthogonal? " << M1.isOrthogonal() << '\n' << '\n';

    // Testing multiplication of two Matrix objects
    Matrix M4 = M1 * M3;
    M4.print();

    // Testing arithmetic operations between Matrix and constant
    Matrix M5 = M1 + 4;
    M5.print();
    Matrix M6 = M1 - 1;
    M6.print();
    Matrix M7 = M1 * 5;
    M7.print();
    Matrix M8 = M1 / 2;
    M8.print();

    // Testing addition and subtraction between two Matrix objects
    Matrix M9 = M1 + M1;
    M9.print();
    Matrix M10 = M1 - M1;
    M10.print();

    return 0;
}
