#include "MethodesCreux.h"
#include "PreconditionneurCreux.h"
#include "LectureCreux.h"
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
	cout << "6) truc a vin100" << endl;
	cin >> choix_matrice;
	system ("clear");

	int kmax = 1000;
	SparseMatrix<double> An;
	std::string source;

	switch(choix_matrice)
	{
		case 1:
		{
			int n;
			kmax = 1000;
			system ("clear");

			cout << "------------------------------------" << endl;
			cout << "Paramètres de la matrice A" << endl;
			cout << "------------------------------------" << endl;
			cout << "n = ?" << endl;
			cin >> n;
			system ("clear");
			An.resize(n,n);

			vector<Triplet<double>> triplets;
			for (int i=0; i<An.rows(); ++i)
			{
				triplets.push_back({i,i,16.});
				if (i > 0)
				triplets.push_back({i,i-1,-2.});
				if (i < An.rows()-1)
				triplets.push_back({i,i+1,-2.});
			}
			An.setFromTriplets(triplets.begin(), triplets.end());

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

		case 6:
			source = "Vincente.mtx";
			Read(An,source);
			break;


		default:
			cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
			exit(0);
	}

	int choix;
	cout << "------------------------------------" << endl;
	cout << "Choississez la méthode : " << endl;
	cout << "1)  Jacobi"<< endl;
	cout << "2)  Gradient à pas Optimal" << endl;
	cout << "3)  Résidu Minimum" << endl;
	cout << "4)  GMRes" << endl;
	cout << "5)  Les 4 en même temps" << endl;
	cout << endl;
	cout << "6)  Résidu minimal préconditionné à gauche" << endl;
	cout << "7)  Résidu minimal préconditionné SSOR" << endl;
	cout << "8)  Résidu minimal préconditionné à droite" << endl;
	cout << "9)  Résidu minimal préconditionné à droite & flexible avec SSOR" << endl;
	cout << "10) Résidu minimal préconditionné à droite & flexible avec le Résidu Minimum" << endl;

	cin >> choix;
	system ("clear");

	double erreur = 1.0e-10;

	SparseVector<double> b(An.rows());
	b = An.col(0);
//  VectorXd c = VectorXd::Random(An.rows()).array().abs();
//	b = c.sparseView();
	SparseVector<double> x0(An.rows());
	x0.setZero();

	if (choix_matrice==1)
	{
		cout << "------------------------------------" << endl;
		cout << "Affichage de la matrice A" << endl;
		cout << "------------------------------------" << endl;
	//	cout << "A = \n" << An << endl;

		cout << "------------------------------------" << endl;
		cout << "Affichage du vecteur b" << endl;
		cout << "------------------------------------" << endl;
	//	cout << "b = \n" << b << endl;

	}

	if (choix_matrice==6)
	{
		b.coeffRef(0) = 0.;
		b.coeffRef(1) = 1.;
		b.coeffRef(2) = 2.;
	}

	cout << "\n" << "------------------------------------" << endl;

	//Définition d'un poiteur de Methodes
	Methodes* sys(0);
	MethodesPreconditionnees* pop(0);
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

		case 6:
			pop =	new RMPJ(An, b, x0, erreur, kmax);
			break;

		case 7:
			pop =	new RMPSSOR(An, b, x0, erreur, kmax);
			break;

		case 8:
			pop =	new RMPD(An, b, x0, erreur, kmax);
			break;

		case 9:
			pop =	new RMPDflexSSOR(An, b, x0, erreur, kmax);
			break;

		case 10:
			pop =	new RMPDflexResMin(An, b, x0, erreur, kmax);
			break;

		default:
			cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
			exit(0);
	}

	if (choix < 6)
	{
		sys->AlgoGene();
	}
	else
	{
		pop->AlgoGene();
	}
	cout << "\n------------------------------------ \n" << endl;


	return 0;
}
