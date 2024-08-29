#include "potential.cpp"


class Point
{
    public:
    double fg; // специально не храню r, cos_th и тд
    double wg;

    Point()
    {
        double r_xe_O2 = sample_cubic(r_xe_O2_min, r_xe_O2_max);
        double cos_th_max;
        if (r_xe_O2 >= 3.857452) // проверено экспериментально
            cos_th_max = 1;
        else
            cos_th_max = dihotomy_solve(0, 1, r_xe_O2);
        double cos_th = sample_linear(-cos_th_max, cos_th_max);
        double energy = -V_xe_O2(r_xe_O2, cos_th);
        if (energy < 0)
            cout << V_xe_O2(r_xe_O2, cos_th) << " " << r_xe_O2 << " " << cos_th << " " << cos_th_max << endl;
        double momentum_O2_max = pow(2 * mu_O2 * energy, 0.5);
        double momentum_O2 = sample_quadratic(0, momentum_O2_max);
        double momentum_xe_O2_max = pow(2 * mu_xe_O2 * (energy - pow(momentum_O2, 2) / (2 * mu_O2)), 0.5);
        double momentum_xe_O2 = sample_cubic(0, momentum_xe_O2_max);
        wg = calculate_weight(momentum_O2_max, momentum_xe_O2_max, cos_th_max);
        fg = exp(-hamiltonian(momentum_O2, momentum_xe_O2, r_xe_O2, cos_th) / T);

    }
    private:
    static double calculate_weight(const double& momentum_O2_max, const double& momentum_xe_O2_max, const double& cos_th_max)
    {
        return 4 * pi * pow(r_O2_e, 2) *
                2 * pi / 3 * (pow(r_xe_O2_max, 3) - pow(r_xe_O2_min, 3)) * (2 * cos_th_max) *
                pi * pow(momentum_O2_max, 2) *
                4. / 3 * pi * pow(momentum_xe_O2_max, 3);
    }
    static double hamiltonian(const double& momentum_O2, const double& momentum_xe_O2, const double& r_xe_O2, const double& cos_th)
    {
        return pow(momentum_O2, 2) / (2 * mu_O2) + pow(momentum_xe_O2, 2) / (2 * mu_xe_O2) + V_xe_O2(r_xe_O2, cos_th);
    }
    static double sample_linear(const double& x_min, const double& x_max)
    {
        return x_min + randomize() * (x_max - x_min);
    }
    static double sample_quadratic(const double& x_min, const double& x_max)
    {
        return pow(pow(x_min, 2) + randomize() * (pow(x_max, 2) - pow(x_min, 2)), 1./2);
    }
    static double sample_cubic(const double& x_min, const double& x_max)
    {
        return pow(pow(x_min, 3) + randomize() * (pow(x_max, 3) - pow(x_min, 3)), 1./3);
    }
    static double randomize()
    {
        return distrib(gen);
    }
    double dihotomy_solve(double a, double b, double r, double error = 1E-8){
        double c = (b + a) / 2;
        int count = 1;
        double V = V_xe_O2(r, c);
        while (b - a > error || V > 0){
            if (V * V_xe_O2(r, a) < 0) b = c; else if (V * V_xe_O2(r, a) > 0) a = c; else {break;}
            c = (b + a) / 2;
            V = V_xe_O2(r, c);
            count++;
        }
    return c;
}

};
class Calculation
{
    unique_ptr<Point[]> points;
    unsigned long long int n;
    public:
    double partition_func;
    double partition_func_error;
    Calculation(unsigned long long int n)
    {
        partition_func = 0;
        partition_func_error = 0;
        this->n = n;
        points = make_unique<Point[]>(n);
    }

    void sum_points()
    {
        for (int i = 0; i < this->n; i++)
            this->partition_func += points[i].fg * points[i].wg;
        this->partition_func = this->partition_func / n;
        for (int i = 0; i < this->n; i++)
            this->partition_func_error += pow(points[i].fg * points[i].wg - this->partition_func, 2);
        this->partition_func = this->partition_func / sym_number / pow(h, 5) * pow((1.408456 * 1E-16),5);
        this->partition_func_error = pow(this->partition_func_error / (this->n * (this-> n - 1)), 0.5) / sym_number / pow(h, 5) * pow((1.408456 * 1E-16),5);
    }
    
};
