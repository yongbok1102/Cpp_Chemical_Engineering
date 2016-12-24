#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace std;

//function to calculate deviation
double calDev(double A, double B, double x[], double y[], int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += pow(y[i] - A*x[i]/(x[i]+ B), 2);
	return sum;
}

//function to calculate minimum deviation
double minDev(double dev[], int n)
{
	double min = dev[0];
	for (int i = 0; i < n; i++)
	{
		if (min > dev[i])
			min = dev[i];
	}
	return min;
}

//function to calculate minimum and maximum
double min(double x[], int n)
{
	double s = x[0];
	for (int i = 0; i < n; i++)
	{
		if (s > x[i])
			s = x[i];
	}
	return s;
}
double max(double x[], int n)
{
	double s = x[0];
	for (int i = 0; i < n; i++)
	{
		if (s < x[i])
			s = x[i];
	}
	return s;
}

int main()
{

//data input
	int n;
	cout << "Enter the number of data\n"; cin >> n;
	double* x; x = new double[n];
	double* y; y = new double[n];
	cout << "Enter the data x_i and y_i\n";
	for (int i = 0; i < n; i++)
	{
		cin >> x[i]; cin >> y[i];
	}

//iteration number
  int itr = 10000000;
  
//y=A*x/(x+B),B = -(y2 - y1) / (y2/x2 - y1/x1), A = y1 + B*y1/x1;
	double A, B;
  
//arbitary two points
	double x1, y1, x2, y2;

	double* setA; setA = new double[itr];
	double* setB; setB = new double[itr];
	double* setDev; setDev = new double[itr];

	int N = 0;
	ofstream out;

	out.open("data.dat");
	for (int i = 0; i < n; i++)
	{
		out << x[i] << '\t' << y[i] << endl;
	}
	out.close();

	srand((unsigned int)time(NULL));
	while (N != itr)
	{
  //generating two points
		x1 = (max(x, n) - min(x, n))*(double)rand() / RAND_MAX + min(x, n);
		y1 = (max(y, n) - min(y, n))*(double)rand() / RAND_MAX + min(y, n);
		x2 = (max(x, n) - min(x, n))*(double)rand() / RAND_MAX + min(x, n);
		y2 = (max(y, n) - min(y, n))*(double)rand() / RAND_MAX + min(y, n);

  //generating the curve passing through the two different points and calculating deviation at the Nth iteration
		if (x1 != x2)
		{
			B = -(y2 - y1) / (y2/x2 - y1/x1);
			A = y1 + B*y1/x1;
			setA[N] = A; setB[N] = B;
			setDev[N] = calDev(A, B, x, y, n);
			N++;
			if (N % 10000 == 0)
				cout << "Itr: " << N << ", Dev: " << setDev[N-1] << endl;
		}
	}

//calculating the minimum deviation
	double devMin = minDev(setDev, itr);

//estimating the line that minimizes the deviation
	out.open("curvefit.txt");
	for (int i = 0; i < itr; i++)
	{
		if (setDev[i] == devMin)
		{
			cout << "y=" << setA[i] << "*x/("<<setB[i]<<"+x), Dev: "<<setDev[i]<<endl;
			out << "y=" << setA[i] << "*x/("<<setB[i]<<"+x), Dev: "<<setDev[i]<<endl;
			break;
		}
	}
	out.close();

	delete[] x;
	delete[] y;
	delete[] setA;
	delete[] setB;
	delete[] setDev;

	return 0;
}
