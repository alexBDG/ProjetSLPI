#include "Methodes.h"
#include "Lecture.h"
#include <string>
#include <iostream>

using namespace std;
using namespace Eigen;

int main()
{
	int choix_matrice;
	system ("clear");
	cout << "------------------------------------" << endl;
	cout << "Quelle matrice utiliser ?" << endl;
	cout << "------------------------------------" << endl;
	cout << "1) Je la crée" << endl;
	cout << "2) BCSSTK18" << endl;
	cout << "3) FIDAPM37" << endl;
	cout << "4) FS 541 4" << endl;
	cout << "5) FS 760 3" << endl;
	cin >> choix_matrice;
	system ("clear");

	int kmax = 50;
	MatrixXd An;
	std::string source;

	switch(choix_matrice)
	{
		case 1:
		{
			int n;
			kmax = 1000;
			double alpha;
			system ("clear");

			cout << "------------------------------------" << endl;
			cout << "Paramètres de la matrice A" << endl;
			cout << "------------------------------------" << endl;
			cout << "n = ?" << endl;
			cin >> n;
			system ("clear");

			alpha = n*n;

			MatrixXd Bn = MatrixXd::Random(n,n).array().abs();
			MatrixXd Bnt = Bn.transpose();
			MatrixXd In = MatrixXd::Identity(n,n);
			An = alpha*In + Bnt*Bn;

			break;
		}

		case 2:
		{
			source = "BCSSTK18/bcsstk18.mtx";
			Read(An,source);
			break;
		}

		case 3:
		{
			source = "FIDAPM37/fidapm37.mtx";
			Read(An,source);
			break;
		}

		case 4:
		{
			source = "FS_541_4/fs_541_4.mtx";
			Read(An,source);
			break;
		}

		case 5:
		{
			source = "FS_760_3/fs_760_3.mtx";
			Read(An,source);
			break;
		}


		default:
			cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
			exit(0);
	}

	int choix;
	cout << "------------------------------------" << endl;
	cout << "Choississez la méthode : " << endl;
	cout << "1) Jacobi"<< endl;
	cout << "2) Gradient à pas Optimal" << endl;
	cout << "3) Résidu Minimum" << endl;
	cout << "4) GMRes" << endl;
	cout << "5) Les 4 en même temps" << endl;
	cin >> choix;
	system ("clear");

	double erreur = 1.0e-10;

	VectorXd b = VectorXd::Random(An.rows()).array().abs();
	VectorXd x0 = VectorXd::Zero(An.rows());

	if (choix_matrice==1)
	{
		cout << "------------------------------------" << endl;
		cout << "Affichage de la matrice A" << endl;
		cout << "------------------------------------" << endl;
		cout << "A = \n" << An << endl;

		cout << "------------------------------------" << endl;
		cout << "Affichage du vecteur b" << endl;
		cout << "------------------------------------" << endl;
		cout << "b = \n" << b << endl;
	}

	cout << "\n" << "------------------------------------" << endl;

	//Définition d'un poiteur de Methodes
	Methodes* sys(0);
	cout << "\n------------------------------------ \n" << endl;


	switch(choix)
	{
		case 1:
			sys = new Jacobi(An, b, x0, erreur, kmax);
			break;

		case 2:
			sys = new GradientOptimal(An, b, x0, erreur, kmax);
			break;

		case 3:
			sys = new ResiduMinimum(An, b, x0, erreur, kmax);
			break;

		case 4:
			sys = new GMRes(An, b, x0, erreur, kmax);
			break;

		case 5:
			sys = new Jacobi(An, b, x0, erreur, kmax);
			sys->AlgoGene();
			sys = new GradientOptimal(An, b, x0, erreur, kmax);
			sys->AlgoGene();
			sys = new ResiduMinimum(An, b, x0, erreur, kmax);
			sys->AlgoGene();
			sys = new GMRes(An, b, x0, erreur, kmax);
			break;

		default:
			cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
			exit(0);
	}

	sys->AlgoGene();
	cout << "\n------------------------------------ \n" << endl;


	return 0;
}
