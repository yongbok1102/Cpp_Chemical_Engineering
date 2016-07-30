#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

/*code for calculating fluid velocity in the pipe */
/*physical properties of fluid */
double rho=1000; double mu=1.005*pow(10.0,-3);
/*characteristic of pipe*/
double L=2.5; double D=0.3; double e=0.26*pow(10.0,-3);
double E=e/D;

double Ploss=8000;

int main()
{
    double u=1;
    double uold=0;
    double tol=pow(10,-6); int N=0;

    ofstream out;
    out.open("pipe.txt");

    double Re,A,B,C;
    double f;

    while(abs(u-uold)>tol){
        uold=u; N+=1;
        Re=uold*rho*D/mu;
        if(Re<2300){
	    /* laminar flow */
            f=64/Re; u=sqrt((2/f)*(D/L)*(Ploss/rho));
            cout<<"iteration "<<N<<" u="<<u<<endl;
            out<<"iteration "<<N<<" u="<<u<<endl;
        }
        else{
	    /*transition and turbulent flow (Serghides's Equation) */
            A=-2*log10(E/3.7 + 12/Re);
            B=-2*log10(E/3.7 + 2.51*A/Re);
            C=-2*log10(E/3.7 + 2.51*B/Re);
            f=pow(A - ((B-A)*(B-A)/(C-2*B+A)),-2);u=sqrt((2/f)*(D/L)*(Ploss/rho));
            cout<<"iteration "<<N<<" u="<<u<<endl;
            out<<"iteration "<<N<<" u="<<u<<endl;
        }
    }

    out.close();
    return 0;
}
