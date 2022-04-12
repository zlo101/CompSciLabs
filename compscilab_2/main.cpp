#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

float mean(float const psi[], float const pdf[], float const dv, unsigned size) {
    float sum = 0.f;
    for (int i = 0; i < size; ++i)
        sum += psi[i] * pdf[i];
    return dv * sum;
}

float mean_recursive(float const psi[], float const pdf[], float const dv, int start, int end) {
    if (start == end)
        return (psi[start] * pdf[start] * dv);
    else if (abs(start - end) == 1)
        return ((psi[start] * pdf[start] + psi[end] * pdf[end]) * dv);
    else {
        auto mid = (int) ((start + end) / 2) ;
        return mean_recursive(psi, pdf, dv, start, mid) + mean_recursive(psi, pdf, dv, mid + 1, end);
    }
}

float mean_close_values(float const psi[], float const pdf[], float const dv, unsigned size) {
    auto* add = new float[size];
    for (int i = 0; i < size; ++i)
        add[i] = psi[i] * pdf[i] * dv;
    int incr = 1;
    while (incr < size) {
        for (int i = 0; i < size - incr; i += 2 * incr)
            add[i] += add[i + incr];
        incr *= 2;
    }
    return add[0];
}

float mean_kahan(float const psi[], float const pdf[], float const dv, unsigned size) {
    float sum = 0.f;
    float t = 0.f;
    for (int i = 0; i < size; ++i) {
        float y = psi[i] * pdf[i] - t;
        float z = sum + y;
        t = (z - sum) - y;
        sum = z;
    }
    return sum * dv;
}

float mean_fma(float const psi[], float const pdf[], float const dv, unsigned size) {
    float sum = 0.f;
    for (int i = 0; i < size; ++i)
        sum = fma(psi[i], pdf[i], sum);
    return dv * sum;
}

double mean_precise(const float psi[], const float pdf[], float const dv, unsigned size) {
    double sum = 0.f;
    for (int i = 0; i != size; ++i)
        sum += (double) psi[i] * pdf[i];
    return (double) dv * sum;
}

void MaxwellDistributionTest(double T) {
    const float dv = 1e-3;
    const int size = 1e6;
    auto* psi = new float[size];
    auto* pdf = new float[size];
    auto* add_psi = new float[size];
    double v = -((double)size / 2) * dv;

    for (int i = 0; i < size; ++i) {
        psi[i] = (float) v;
        add_psi[i] = (float) abs(v);
        pdf[i] = 1.0 / sqrt(T * M_PI) * exp(-(float)v * (float)v / T);
        v += dv;
    }

    cout << "Maxwell distribution for speed values: \n";
    cout << "Normal:                  " << fixed << setprecision(10) << mean(psi, pdf, dv, size) << '\n';
    cout << "Recursive:               " << fixed << setprecision(10) <<  mean_recursive(psi, pdf, dv, 0, size - 1 ) << '\n';
    cout << "Close value sums:        " << fixed << setprecision(10) << mean_close_values(psi, pdf, dv, size) << '\n';
    cout << "Kahan sums:              " << fixed << setprecision(10) << mean_kahan(psi, pdf, dv, size) << '\n';
    cout << "FMA sums:                " << fixed << setprecision(10) << mean_fma(psi, pdf, dv, size) << '\n';
    cout << "Precise value sums:      " << fixed << setprecision(10) << mean_precise(psi, pdf, dv, size) << '\n';
    cout << '\n';
    cout << "Maxwell distribution for absolute speed values: \n";
    cout << "Theoretical prediction:  " << fixed << setprecision(10) << sqrt(T / M_PI) << '\n';
    cout << "Normal:                  " << fixed << setprecision(10) << mean(add_psi, pdf, dv, size) << '\n';
    cout << "Recursive:               " << fixed << setprecision(10) <<  mean_recursive(add_psi, pdf, dv, 0, size - 1 ) << '\n';
    cout << "Close value sums:        " << fixed << setprecision(10) << mean_close_values(add_psi, pdf, dv, size) << '\n';
    cout << "Kahan sums:              " << fixed << setprecision(10) << mean_kahan(add_psi, pdf, dv, size) << '\n';
    cout << "FMA sums:                " << fixed << setprecision(10) << mean_fma(add_psi, pdf, dv, size) << '\n';
    cout << "Precise value sums:      " << fixed << setprecision(10) << mean_precise(add_psi, pdf, dv, size) << '\n';

}

int main() {
    double T;
    cout << "Type in the value of T: ";
    cin >> T;
    MaxwellDistributionTest(T);
}