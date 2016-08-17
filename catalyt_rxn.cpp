#include <iostream>
#include <cmath>
#include <fstream>

/*-------------------------------------------------*/
/*Source code for catalytic reaction	           */
/*Written by yongbok1102			   */
/*-------------------------------------------------*/

using namespace std;

//residual
double Res(double X[], double Xold[], int N)
{
	double sum = 0;
	for (int i = 0; i<N; i++)
	{
		sum += abs(X[i] - Xold[i]);
	}
	return sum / N;
}

//midpoint
double mid(double x, double y)
{
	return (x + y) / 2;
}

int main()
{
	//parameters
	double phi = 1.1; double gamma = 30; double beta = 0.12;

	//mesh generation
	int N = 25; double dr = 1. / N;
	double* R; R = new double[N + 1];
	for (int i = 0; i <= N; i++)
	{
		R[i] = i*dr;
	}

	//initialization
	double* C; C = new double[N + 1];
	double* T; T = new double[N + 1];
	for (int i = 0; i <= N; i++)
	{
		C[i] = 1;
		T[i] = 1;
	}
	double* Cold; Cold = new double[N + 1];
	double* Told; Told = new double[N + 1];

	double dt = 0.00001;
	double resC = 1; double resT = 1;
	int itr = 0;
	ofstream out;
	
	out.open("residual.dat");
	out<<"itr"<<'\t'<< "resC" <<'\t'<<"resT"<<endl;
	while (resC > 1e-013 && resT > 1e-013)
	{
		//initialization
		for (int i = 0; i <= N; i++)
		{
			Cold[i] = C[i];
			Told[i] = T[i];
		}

		//Euler method
		for (int i = 1; i < N; i++)
		{
			C[i] = Cold[i]
				+ dt*pow(dr*R[i], -2)*(pow(mid(R[i], R[i + 1]), 2)*Cold[i + 1] - (pow(mid(R[i], R[i + 1]), 2) + pow(mid(R[i], R[i - 1]), 2))*Cold[i] + pow(mid(R[i], R[i - 1]), 2)*Cold[i - 1])
				- dt*pow(phi, 2)*Cold[i] * exp(gamma*(1 - pow(Told[i], -1)));
			T[i] = Told[i]
				+ dt*pow(dr*R[i], -2)*(pow(mid(R[i], R[i + 1]), 2)*Told[i + 1] - (pow(mid(R[i], R[i + 1]), 2) + pow(mid(R[i], R[i - 1]), 2))*Told[i] + pow(mid(R[i], R[i - 1]), 2)*Told[i - 1])
				+ dt*beta*pow(phi, 2)*Cold[i] * exp(gamma*(1 - pow(Told[i], -1)));
		}
		
		//boundary condition
		C[N]=1; T[N]=1;
		C[0]=C[1]; T[0]=T[1];
		itr++;
		resC = Res(C, Cold, N + 1); resT = Res(T, Told, N + 1);
		cout << "iteration: " << itr << " residual(C): " << resC << " residual(T): " << resT << endl;
		out<<itr<<'\t'<<resC<<'\t'<<resT<<endl;
		if(itr==5000000)
			break;
	}
	out.close();
	out.open("proj0007_concentration.dat");
	for(int i=0;i<=N;i++)
	{
		out<<R[i]<<'\t'<<C[i]<<endl;
	}
	out.close();
	
	out.open("proj0007_temperature.dat");
	for(int i=0;i<=N;i++)
	{
		out<<R[i]<<'\t'<<T[i]<<endl;
	}
	out.close();
	
	delete[] R;
	delete[] C; delete[] Cold;
	delete[] T; delete[] Told;

	return 0;
}
