double min(double x[], int N)
{
    double min=x[0];
    for(int i=0;i<N;i++)
    {
        if(min<x[i])
            min=min;
        else
            min=x[i];
    }
    return min;
}

double max(double x[],int N)
{
    double max=x[0];
    for(int i=0;i<N;i++)
    {
        if(max>x[i])
            max=max;
        else
            max=x[i];
    }
    return max;
}
