#include <fstream>
#include <string>
#include "LectureCreux.h"

using namespace Eigen;
using namespace std;

void Read(SparseMatrix<double> & matrix, string & source)
{
  //Create the file:
  std::string file;
  file = "Matrices/" + source;

  // Open the file:
  std::ifstream fin(file);
  // Declare variables:
  int M, N, L;

  // Ignore headers and comments:
  while (fin.peek() == '%') fin.ignore(2048, '\n');

  // Read defining parameters:
  fin >> M >> N >> L;

  // Create your matrix:
  matrix.resize(M,N);
  matrix.setZero();	     // Creates the array of M*N size

  // Read the data
  for (int l = 0; l < L; l++)
  {
    int m, n;
    double data;
    fin >> m >> n >> data;
    matrix.coeffRef(m-1,n-1) = data;
    //matrix.coeffRef(n-1,m-1) = data;
  }

  fin.close();
}
