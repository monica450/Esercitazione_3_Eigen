#include <iostream>
#include "Eigen/Eigen"
#include <iomanip>

using namespace std;
using namespace Eigen;

int main()
{   
    Vector2d x_esatta;
	x_esatta<<-1.0e+0,-1.0e+00;
    cout << scientific << setprecision(1);
    Matrix2d A1;
	A1 <<5.47001962252291e-01,-3.770900990025203e-02, 8.320502943378437e-01,-9.992887623566787e-01;
	
	Vector2d b1;
	b1<<-5.169911863249772e-01, 1.672384680188350e-01;
	
	//PALU
	PartialPivLu<Matrix2d> palu1(A1); //esegue la fattorizzazione
	Vector2d x1_palu=palu1.solve(b1);  //risolve il sistema col vettore b
	
	//QR
	HouseholderQR<Matrix2d> qr1(A1);
	Vector2d x1_qr=qr1.solve(b1);
	
	//SOLUZIONE
	cout <<"Soluzione PALU:" << x1_palu << endl;
	cout <<"Soluzione QR:" << x1_qr << endl;
	
	double err_rel_palu1;
    err_rel_palu1<<	(x1_palu-x_esatta).norm()/(x_esatta.norm());
	cout<<err_rel_palu1<<endl;


    Matrix2d A2;
	A2 <<5.547001962252291e-01,-5.540607316466765e-01, 8.320502943378437e-01,-8.324762492991313e-01;
	
	Vector2d b2;
	b2<<-6.394645785530173e-04, 4.259549612877223e-04;
	
	//PALU
	PartialPivLu<Matrix2d> palu2(A2); //esegue la fattorizzazione
	Vector2d x2_palu=palu2.solve(b2);  //risolve il sistema col vettore b
	
	//QR
	HouseholderQR<Matrix2d> qr2(A2);
	Vector2d x2_qr=qr2.solve(b2);
	
	//SOLUZIONE
	cout <<"Soluzione PALU:" << x2_palu << endl;
	cout <<"Soluzione QR:" << x2_qr << endl;

    double err_rel_palu2;
    err_rel_palu2<<	(x2_palu-x_esatta).norm()/(x_esatta.norm());
	cout<<err_rel_palu2<<endl;
	
	Matrix2d A3;
	A3 <<5.547001962252291e-01,-5.547001955851905e-01, 8.320502943378437e-01,-8.320502947645361e-01;
	
	Vector2d b3;
	b3<<-6.400391328043042e-10, 4.266924591433963e-10;
	
	//PALU
	PartialPivLu<Matrix2d> palu3(A3); //esegue la fattorizzazione
	Vector2d x3_palu=palu3.solve(b3);  //risolve il sistema col vettore b
	
	//QR
	HouseholderQR<Matrix2d> qr3(A3);
	Vector2d x3_qr=qr3.solve(b3);
	
	//SOLUZIONE
	cout <<"Soluzione PALU:" << x3_palu << endl;
	cout <<"Soluzione QR:" << x3_qr << endl;

    Vector2d err_rel_palu3;
    err_rel_palu3<<	(x3_palu-x_esatta).norm()/(x_esatta.norm());
	cout<<err_rel_palu3<<endl;
	

    return 0;
}
