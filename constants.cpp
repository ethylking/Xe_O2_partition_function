#include "libraries.cpp"

#define q_el 1.6E-12 //  ESU 
#define T (300 * 1.380649 * 1E-16 / q_el * 8065.54) // cm^-1 
#define h (6.62607015 * 1E-27) // plank constant
#define pi 3.14159265359

// Global constants
double r_xe_O2_min = 3.4352301326015; //  Angstrom 
double r_xe_O2_max = 12.; //  Angstrom

static const double m_xe = 131.293 * 1.6605402 * 1E-24; // g
static const double m_O2 = 32 * 1.6605402 * 1E-24;

static const double r_O2_e = 1.20752; //    Angstrom
static const double I_O2 = m_O2 / 2 * pow((r_O2_e * 1E-8), 2);
static const double mu_O2 = m_O2 / 2;
static const double mu_xe_O2 = m_O2 * m_xe / (m_xe + m_O2);

static const double v_O2 = 1580.; //  O2 vibration quantum, cm^-1
static const int sym_number = 2; // symmetry caused by oxygen
static const double Z_xe = ((m_xe * T / 8065.54 * q_el / pow(h, 2)), 1.5);
static const double Z_O2 = ((m_O2 * T / 8065.54 * q_el / pow(h, 2)), 1.5) * 
        (T / sym_number / pow(h / (2 * pi), 2) * 2 * I_O2 / 8065.54 * q_el);
static const double Z_xe_O2 = 0;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<double> distrib(0, 1);
