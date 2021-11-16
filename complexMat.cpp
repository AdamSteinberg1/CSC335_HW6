//Problem 2 from Dr. Pounds

#include <iostream>
#include <cmath>
#include <armadillo>

using namespace arma;
using namespace std;

//return a nxn matrix M
cx_mat make_M(int n)
{
  cx_mat M(n,n);
  for(int j = 1; j <= n; j++)
  {
    for (int k = 1; k <= n; k++)
    {
      double term1 = cos(j*j*k*M_PI/pow(n,3));
      double term2 = sin(2*j*j*k*M_PI/pow(n,3));
      double alpha = pow(pow(term1, 2) + pow(term2, 2), -0.5);
      complex<double> entry(alpha*term1, alpha*term2);
      M(j-1,k-1) = entry;
    }
  }
  return M;
}

//sums the elements along a matrix's diagonal
complex<double> diagonalSum(cx_mat M)
{
  complex<double> sum = 0.0;
  for(auto& entry : cx_vec(M.diag()))
  {
    sum += entry;
  }
  return sum;
}

int main()
{
  //part a
  //generate a table of dimension and diagonal sum
  cout << "N\tDiagonal Sum" << endl;
  for(int i = 1; i<=100; i++)
  {
    cx_mat M = make_M(i);
    cx_mat P = M*M.i();
    cout << i << '\t' << diagonalSum(P) << endl;
  }

  //part b find the largest dimension for which the multiplication is possible

  //Armadillo can multiply matrices of any size, provided there is sufficient
  //memory. This means the largest dimension for which the multiplication is
  //possible, depends on the current load on Cobra. I've found that usually
  //there is enough memory available to multiply 15000x15000 matrices. Much
  //larger and I get an out of memory error.

  //this will succeed if Cobra is not under heavy load
  const int MAX_DIM = 20000;
  cx_mat M1 = make_M(MAX_DIM);
  cx_mat P1 = M1*M1.i();
  cout << MAX_DIM << '\t' << diagonalSum(P1) << endl;


  //this will fail with the following error
  //error: arma::memory::acquire(): out of memory
  //terminate called after throwing an instance of 'std::bad_alloc'
  //  what():  std::bad_alloc

  cx_mat M2 = make_M(50000);
  cx_mat P2 = M2*M2.i();
  cout << 50000 << '\t' << diagonalSum(P2) << endl;

}
