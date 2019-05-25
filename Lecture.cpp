#include <fstream>
#include "Lecture.h"

using namespace Eigen;
using namespace std;

void Read(MatrixXd & matrix, string & source)
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
  matrix = MatrixXd::Zero(M,N);	     // Creates the array of M*N size

  // Read the data
  for (int l = 0; l < L; l++)
  {
  int m, n;
  double data;
  fin >> m >> n >> data;
  matrix(m-1,n-1) = data;
  }

  fin.close();
}
