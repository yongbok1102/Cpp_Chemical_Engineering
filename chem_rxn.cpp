#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double mid(double x, double y)
{
	return (x+y)/2;
}

int main()
{
	//note : i: z-dir, j:r-dir
	
	//parameters for this problem
	double Pe = 0.5; double Da = 10;
	
	//grid generation
	int nr = 10; int nz = 1000;
	int np = (nr+1)*(nz+1);
	int ncell = nr*nz;
	ofstream out;
	out.open("proj0004.vtk");
	out << "# vtk DataFile Version 3.1\n";
	out << "chemical reaction\n";
	out << "ASCII\n";
	out << "DATASET UNSTRUCTURED_GRID\n";
	out << "POINTS " << np << " float\n";	
	double* R; double* Z;
	R = new double[np]; Z = new double[np];
	double dr = 1./nr; double dz = 1./nz;
	int n=0;
	for(int i=0;i<=nz;i++)
	{
		for(int j=0;j<=nr;j++)
		{
			R[n] = j*dr;
			Z[n] = i*dz;
			n++;
		}
	}
	for(int i=0;i<np;i++)
	{
		out<<Z[i]<<'\t'<<R[i]<<'\t'<<0<<endl;
	}
	out << "CELLS " << ncell << '\t' << ncell * 5 << endl;
	for(int i=0;i<np;i++)
	{
		if(i%(nz+1)!=nz && i/(nz+1)!=nr)
			out<<4<<'\t'<<i<<'\t'<<i+1<<'\t'<<i+1+(nz+1)<<'\t'<<i+nz+1<<endl;
	}
	out << "CELL_TYPES " << np << endl;
	for (int i = 0; i < np; i++)
	{
		out << 9 << endl;
	}

	double* F; double* Fold;
	F = new double[nr + 1]; Fold = new double[nr + 1];

	//Inlet value
	for (int i = 0; i <= nr; i++)
	{
		F[i] = 1;
	}
	out << "POINT_DATA " << np << endl;
	out << "SCALARS F float\n";
	out << "LOOKUP_TABLE default\n";
	for (int i = 0; i <= nr; i++)
	{
		out << F[i] << endl;
	}
	for (int j = 1; j<= nz; j++)
	{
		//Required for the forward difference schemes of dF/dz
		for (int i = 0; i <= nr; i++)
		{
			Fold[i] = F[i];
		}
		for (int i = 0; i <= nr; i++)
		{
			if (i != 0 && i != nr)
				F[i] = Fold[i] + dz*(pow(R[i] * dr*dr*Pe, -1)*(mid(R[i], R[i + 1])*(Fold[i + 1] - Fold[i]) - mid(R[i], R[i - 1])*(Fold[i] - Fold[i - 1])) - Da*Fold[i]);
			//boundary condition
			else if (i == nr)
				F[i] = F[i - 1];
			else
				F[i] = F[i + 1];
		}
		for (int i = 0; i <= nr; i++)
		{
			out << F[i] << endl;
		}
		cout << "calculation " << j << endl;
	}
	out.close();
	return 0;
}
