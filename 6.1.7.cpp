//Problem 7 from section 6.1

#include <iostream>
#include <armadillo>
#include <cassert>

using namespace arma;
using namespace std;


//calculates the vector that solves an augmented matrix
//matrix is a n x n+1 augmented matrix
vec gaussianElimination(mat matrix)
{
  int n = matrix.n_rows;
  assert(n+1 == matrix.n_cols);

  for(int i = 0; i < n-1; i++)
  {
    float p;
    for(int k = i; k < n; k++)
    {
      if(matrix(k, i) != 0.0)
      {
        p = k;
        break;
      }
      else if(k==n-1)
      {
        cerr << "no unique solution exists" << endl;
        return vec();
      }
    }

    if(p!=i)
    {
      matrix.swap_rows(p,i);
    }

    for(int j = i+1; j < n; j++)
    {
      float m = matrix(j,i) / matrix(i,i);
      matrix.row(j) = matrix.row(j) - m * matrix.row(i);
    }
  }

  if(matrix(n-1,n-1) == 0)
  {
    cerr << "no unique solution exists" << endl;
    return vec();
  }

  vec x(n);
  x(n-1) = matrix(n-1,n) / matrix(n-1, n-1);
  for(int i = n-2; i >= 0; i--)
  {
    float sum = 0;
    for(int j = i+1; j < n; j++)
    {
      sum += matrix(i,j)*x(j);
    }
    x(i) = (matrix(i,n) - sum) / matrix(i,i);
  }

  return x;
}

int main()
{
  cout.precision(10);
  //uses two methods to solve the system of linear equations
  //first method is a custom function to perform Gaussian Elimination
  //based on algorithm 6.1 from Burden text
  mat A =  {{1.0,     1.0/2.0,  1.0/3.0,  1.0/4.0,  1.0/6.0},
            {1.0/2.0, 1.0/3.0,  1.0/4.0,  1.0/5.0,  1.0/7.0},
            {1.0/3.0, 1.0/4.0,  1.0/5.0,  1.0/6.0,  1.0/8.0},
            {1.0/4.0, 1.0/5.0,  1.0/6.0,  1.0/7.0,  1.0/9.0}};
  vec x = gaussianElimination(A);
  for(int i = 0; i < x.size(); i++)
    cout << "x" << i+1 << " = " << x(i) << endl;
  cout << endl;
  //second method is armadillo's built-in solve function
  A = {{1.0,     1.0/2.0,  1.0/3.0,  1.0/4.0},
       {1.0/2.0, 1.0/3.0,  1.0/4.0,  1.0/5.0},
       {1.0/3.0, 1.0/4.0,  1.0/5.0,  1.0/6.0},
       {1.0/4.0, 1.0/5.0,  1.0/6.0,  1.0/7.0}};

  vec b = {1.0/6.0, 1.0/7.0, 1.0/8.0, 1.0/9.0};
  x = solve(A, b);
  for(int i = 0; i < x.size(); i++)
    cout << "x" << i+1 << " = " << x(i) << endl;

}
