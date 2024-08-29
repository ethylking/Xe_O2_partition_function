#include "classes.cpp"

int main()
{
    unsigned long long int n = 5E5;
    cout << setprecision(16);
    Calculation experiment(n);
    experiment.sum_points();
    cout << "partition func is " << experiment.partition_func 
    << "\t" << "error is " << 
    experiment.partition_func_error << endl <<
    "relative error is " << experiment.partition_func_error / experiment.partition_func << endl;

    return 0;
}