#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "minmax.h"
using namespace std;

double rA(double C)
{
    static double k=0.45;
    static double n=1.5;
    return k*pow(C,n);
}

double f(double C)
{
    static double Cin=100;
    static double Cout=10;
    return (Cin-C)*(1/(rA(Cout))-1/(rA(C)));
}

int main()
{
    double a=10; double b=100;
    int N;
    cout<<"Enter the iteration number\n";cin>>N;
    double x;
    
    double* pF;
    pF = new double[N];
    
    srand((unsigned int)time(NULL));
    for(int i=0;i<N;i++)
    {
        x=(b-a)*(double)rand()/RAND_MAX + a;
        pF[i]=f(x);
        cout<<i+1<<'\t'<<pF[i]<<endl;
    }
    
    double res=max(pF,N);
    
    double tol=1e-010;
    double C = (a+b)/2; double rel;
    cout<<"Enter the relaxation factor\n";cin>>rel;
    while(abs(f(C)-res)>tol){
        C=C-rel*(f(C)-res);
        cout<<"C="<<C<<endl;
    }
    double F0;
    cout<<"Enter the inlet flow rate[L/min]\n";cin>>F0;
    
    double V1 = F0*(b-C)/rA(C);
    double V2 = F0*(C-a)/rA(a);
    cout<<"V1="<<V1<<'\t'<<"V2="<<V2<<endl;
    
    delete[] pF;
    return 0;
}
