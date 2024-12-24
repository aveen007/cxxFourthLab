// cxxFourthLab.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "sparseVector.h"
#include "Matrix.h"
using namespace std;
int main()
{
    SparseVector <double> a(10);
    a.set(1, 2);
    a.set(2, 3);
    //cout << a.get(1) << a.getSize() << a[1];  
    //a.set(1, 0);
    //cout << a.get(1) << a.getSize() << a[1];
    SparseVector <double> b(10);
    b.set(2, 4);
    //cout << b.get(1) << b.getSize() << a[1];
   SparseVector<double>c= a + b;
   c.print();
   int d = a.dot( b);
   cout<<endl << d<<endl;

   SparseMatrix<double > myMatrix(2, 2);
   myMatrix.set(0, 0, 2);
   myMatrix.set(1, 1, 2);
   std::cout << myMatrix(1, 1);
   SparseMatrix<double > invertedMatrix= myMatrix.inverse();
   myMatrix = myMatrix ^ 2;
   myMatrix = myMatrix.exp();
   myMatrix=myMatrix.pow(2.5);

   
    std::cout << "Hello World!\n";
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
