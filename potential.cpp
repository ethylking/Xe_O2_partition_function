#include "constants.cpp"

// Potential constants
static const double r_xe_O2_e = 3.87; //  Angstrom
static const double e_xe_O2 = 15.20 * 8.06554; //  Angstrom
static const double b = 6.5;
static const double C0 = 1.02 * 1E5 * 8.06554;
static const double b1 = -2.5457689451844066;
static const double b2 = 5.665376678222135;
static const double b3 = -4.276960450638895;
static const double b4 = -2.2417153960400404;
static const double A2 = 8.4 * 1E6 * 8.06554;
static const double alpha = 3.37;
static const double a2 = 0.26;
static const double A4 = 1E6 * 8.06554;
static const double rm = 4.05;
static const double x1 = 1.12;
static const double x2 = 1.55;

double V01(const double& r)
{
    double x = r / rm;
    return e_xe_O2 * (exp(- 2 * b * (x - 1)) - 2 * exp(-b * (x - 1)));
}
double V02(const double& r)
{
    double x = r / rm;
    return e_xe_O2 * (b1 + b2 * (x - x1) + (x - x2) * (b3 + (x - x1) * b4));
}
double V03(const double& r)
{
    return -C0 / pow(r, 6);
}
double V2(const double& r)
{
    return A2 * exp(-alpha * r) - a2 * C0 / pow(r, 6);
}
double V4(const double& r)
{
    return A4 * exp(-alpha * r);
}
double P2(const double& x)
{
    return 0.5 * (3. * x * x - 1.);
}
double P4(const double& x)
{
    return 1 / 8. * (35. * pow(x, 4) - 30. * x * x + 3.);
}

double V_xe_O2(const double r = r_xe_O2_e, const double th = 0)
{
    double x = r / rm;
    if (x <= x1){
        return V01(r) + V2(r) * P2(th) + V4(r) * P4(th);
    } else if (x > x2){
        return V03(r) + V2(r) * P2(th) + V4(r) * P4(th);
    }
    return V02(r) + V2(r) * P2(x) + V4(r) * P4(x);
}
