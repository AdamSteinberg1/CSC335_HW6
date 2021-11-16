//Problem 1 from Dr. Pounds

#include <armadillo>
#include <iostream>
#include <cmath>

using namespace arma;
using namespace std;

mat hilbert(int n)
{
  mat H(n,n);
  for (int i = 1; i <= n; i++)
  {
    for (int j = 1; j <= n; j++)
    {
      H(i-1,j-1) = 1.0/(i+j-1.0);
    }
  }
  return H;
}

int main()
{
  cout.precision(15);
  double determinant = 1;
  int i = 0;
  while(!isnan(determinant) && determinant != 0.0)
  {
    i++;
    mat H = hilbert(i);
    determinant = det(H);
    cout << "determinant of the " << i << "x" << i << " Hilbert matrix = " << determinant << endl;
  }

  cout << "\nCould not calculate the determinant of the " << i << "x" << i <<" Hilbert matrix." << endl;
  cout << "Therefore, the largest Hilbert matrix that I can calculate the determinant of is " << i-1 << "x" << i-1 << endl;
}
