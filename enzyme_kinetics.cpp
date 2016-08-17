/*-------------------------------------------------*/
/*Source code for enzyme kinetics					*/
/*Written by yongbok1102							*/
/*-------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

//residual
double Res(double X[], double Xold[], int N)
{
	double sum = 0;
	for (int i = 0; i<N; i++)
	{
		sum += pow(X[i] - Xold[i], 2);
	}
	return sqrt(sum / N);
}

//midpoint
double mid(double x, double y)
{
	return (x + y) / 2;
}

int main()
{
	//parameters
	double De = 1e-006;
	double Rmax = 0.50;
	double Km = 3.0*pow(10, -4);
	double r = 0.003;
	double Cout = 0.5;

	//mesh generation
	int N = 600; double dr = r / N;
	double* R; R = new double[N + 1];
	for (int i = 0; i <= N; i++)
	{
		R[i] = i*dr;
	}

	//initialization
	double* C; C = new double[N + 1];
	for (int i = 0; i <= N; i++)
	{
		C[i] = Cout;
	}
	double* Cold; Cold = new double[N + 1];

	int itr = 0; double res = 1;
	ofstream out;
	out.open("enzyme_kinetics_residual.dat");

	//Calculated iteratively and printing out the residual
	while(log10(res) > -40)
	{
		for (int i = 0; i <= N; i++)
		{
			Cold[i] = C[i];
		}
		
		//Gauss-Seidel iteration
		for (int j = 1; j <= 10000; j++)
		{
			for (int i = 0; i <= N; i++)
			{
				if (i == 0)
					C[i] = C[i + 1];
				else if (i == N)
					C[i] = Cout;
				else
					C[i] = (De*pow(R[i] * dr, -2)*pow(mid(R[i], R[i + 1]), 2)*C[i + 1] + De*pow(R[i] * dr, -2)*pow(mid(R[i], R[i - 1]), 2)*C[i - 1]) /
					(De*pow(R[i] * dr, -2)*(pow(mid(R[i], R[i + 1]),2) + pow(mid(R[i], R[i - 1]),2)) + Rmax / (Km + Cold[i]));
			}
		}
		itr++;
		res = Res(C, Cold, N + 1);
		cout.precision(40);
		out.precision(40);
		cout << "iteration: " << itr << " residual: " << res << endl;
		out << itr << '\t' << res << endl;
		if (itr == 1000000)
			break;
	}
	out.close();

	//printing out the result
	out.open("enzyme_kinetics_result.dat");
	for(int i=0;i<=N;i++)
	{
		out<<R[i]<<'\t'<<C[i]<<endl;
	}
	out.close();
	delete[] R;
	delete[] C;
	delete[] Cold;

	return 0;
}