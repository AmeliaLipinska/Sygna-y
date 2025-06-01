#define _USE_MATH_DEFINES 
#include <matplot/matplot.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "sygnalytest.h"
#include <string>

using namespace matplot;
using namespace std;

std::vector<double> operator*(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());
    std::transform(vec.begin(), vec.end(), result.begin(),
        [scalar](double v) { return v * scalar; });
    return result;
}

// Operator do mnozenia skalar przez wektor (komutatywnosc)
std::vector<double> operator*(double scalar, const std::vector<double>& vec) {
    return vec * scalar; // reuse the other overload
}

// --- Generowanie sygnalow ---
vector<double> generate_sine(double freq, double t_start, double t_end, size_t num_samples) {
    vector<double> x = linspace(t_start, t_end, num_samples);
    vector<double> y(num_samples);
    for (size_t i = 0; i < num_samples; ++i) {
        y[i] = sin(2 * M_PI * freq * x[i]);
    }
    return y;
}

vector<double> generate_cosine(double freq, double t_start, double t_end, size_t num_samples) {
    vector<double> x = linspace(t_start, t_end, num_samples);
    vector<double> y(num_samples);
    for (size_t i = 0; i < num_samples; ++i) {
        y[i] = cos(2 * M_PI * freq * x[i]);
    }
    return y;
}

vector<double> generate_square(double freq, double t_start, double t_end, size_t num_samples) {
    vector<double> x = linspace(t_start, t_end, num_samples);
    vector<double> y(num_samples);
    double period = 1.0 / freq;
    for (size_t i = 0; i < num_samples; ++i) {
        y[i] = (fmod(x[i], period) < period / 2.0) ? 1.0 : -1.0;
    }
    return y;
}

vector<double> generate_sawtooth(double freq, double t_start, double t_end, size_t num_samples) {
    vector<double> y(num_samples);

    if (freq <= 0 || num_samples < 2 || t_end <= t_start) {
        cerr << "Bledne parametry sygnalu!\n";
        return y;
    }

    double T = 1.0 / freq;
    double dt = (t_end - t_start) / (num_samples - 1);

    for (size_t i = 0; i < num_samples; ++i) {
        double t = t_start + i * dt;
        double phase = fmod(t, T) / T;  // faza z przedzialu [0, 1)
        y[i] = 2.0 * phase - 1.0;
    }

    return y;
}


//---------------------------------------------------------------------------------------------------------------------------------------------

//dft odczytuje nadawany impuls i przeksztalca go w postac wektorow, umieszcza odpowiednie wartosci w odpowiednich wektorach
//DFT
DFTResult computeDFT(const vector<double>& signal) {
    size_t N = signal.size();
    DFTResult dft;
    dft.real.resize(N);
    dft.imag.resize(N);

    //dla kazdej czestotliwosci f po kolei sprawdza jak bardzo przypomina sin/cos
    for (size_t f = 0; f < N; ++f) {
        double sum_real = 0.0;
        double sum_imag = 0.0;
        for (size_t n = 0; n < N; ++n) {
            double angle = -2.0 * M_PI * f * n / N;
            sum_real += signal[n] * cos(angle); //zbiera informacje o cosinusie-RE
            sum_imag += signal[n] * sin(angle); //zbiera informacje o sinusie-IM
        }

        //koncowe skladowe dla danej czestotliwosci dla cos i sin
        dft.real[f] = sum_real;
        dft.imag[f] = sum_imag;
    }
    return dft;
}

//idft czyli odwrotna transformata fouriera, przeksztalca wektory z poprzednio odczytanego dft w postac sygnalu
//IDFT
vector<double> computeIDFT(const DFTResult& dft) {
    size_t N = dft.real.size();
    vector<double> signal(N);

    //dla kazdej probki czasowej n
    for (size_t n = 0; n < N; ++n) {
        double sum = 0.0;
        for (size_t k = 0; k < N; ++k) {
            double angle = 2.0 * M_PI * k * n / N;
            sum += dft.real[k] * cos(angle) - dft.imag[k] * sin(angle); //skladanie spowrotem sygnalu z dft
        }

        //nie wiem po chuja to ale chat kaze bo bez tego bedzie rozmazane czy cos(pozniej doczytam)
        signal[n] = sum / N;
    }
    return signal;
}

// --- Wizualizacja DFT (amplituda) ---

void plotDFT(const DFTResult& dft) {
    vector<double> amplitudes(dft.real.size());
    for (size_t i = 0; i < dft.real.size(); ++i) {
        amplitudes[i] = sqrt(dft.real[i] * dft.real[i] + dft.imag[i] * dft.imag[i]);
    }
    plot(amplitudes);
    title("DFT (Amplituda)");
    //show();
}

//Filtracja 1D
vector<double> filter_1D(const vector<double>& signal, const vector<double>& kernel) {
    size_t N = signal.size();
    size_t K = kernel.size();
    size_t half_K = K / 2;

    vector<double> filtered(N, 0.0);

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < K; ++j) {
            int idx = static_cast<int>(i) + static_cast<int>(j) - static_cast<int>(half_K);
            if (idx >= 0 && idx < N) {
                filtered[i] += signal[idx] * kernel[j];
            }
        }
    }
    return filtered;
}
//filtracja 2d
vector<vector<double>> Filter2D(const vector<vector<double>>& signal, const vector<vector<double>>& filter) {
    int filterRows = filter.size();
    int filterCols = filter[0].size();
    int halfRows = filterRows / 2;
    int halfCols = filterCols / 2;

    int rows = signal.size();
    int cols = signal[0].size();
    vector<vector<double>> result(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double sum = 0.0;
            for (int m = -halfRows; m <= halfRows; ++m) {
                for (int n = -halfCols; n <= halfCols; ++n) {
                    int y = i + m;
                    int x = j + n;
                    if (y >= 0 && y < rows && x >= 0 && x < cols) {
                        sum += signal[y][x] * filter[m + halfRows][n + halfCols];
                    }
                }
            }
            result[i][j] = sum;
        }
    }

    return result;
}

// --- Main ---
void plotting(const std::string& type1, double freq1, double amp1,
    const std::string& type2, double freq2, double amp2,
    double t_start, double t_end, size_t num_samples) {

    figure();
    // Generowanie sygnalu
    vector<double> x = linspace(t_start, t_end, num_samples);
    vector<double> y1, y2;

    //sygnal 1 wraz ze zmiana amplitudy
    if (type1 == "sin") y1 = generate_sine(freq1, t_start, t_end, num_samples) * amp1;
    else if (type1 == "cos") y1 = generate_cosine(freq1, t_start, t_end, num_samples) * amp1;
    else if (type1 == "square") y1 = generate_square(freq1, t_start, t_end, num_samples) * amp1;
    else if (type1 == "sawtooth") y1 = generate_sawtooth(freq1, t_start, t_end, num_samples) * amp1;
    else {
        cerr << "Nieznany typ sygnalu 1!\n";
        return;
    }

    //sygnal 2 wraz ze zmiana amplitudy
    if (type2 == "sin") y2 = generate_sine(freq2, t_start, t_end, num_samples) * amp2;
    else if (type2 == "cos") y2 = generate_cosine(freq2, t_start, t_end, num_samples) * amp2;
    else if (type2 == "square") y2 = generate_square(freq2, t_start, t_end, num_samples) * amp2;
    else if (type2 == "sawtooth") y2 = generate_sawtooth(freq2, t_start, t_end, num_samples) * amp2;
    else {
        cerr << "Nieznany typ sygnalu 2!\n";
        return;
    }

    // Sprawdzenie d³ugoœci
    if (y1.size() != y2.size()) {
        cerr << "Sygnaly maj¹ rozne dlugosci!\n";
        return;
    }

    // Sygnal zlozony (suma)
    vector<double> y_combined(y1.size());
    for (size_t i = 0; i < y1.size(); ++i) {
        y_combined[i] = y1[i] + y2[i];
    }

    // DFT
    DFTResult dft = computeDFT(y_combined);

    // IDFT
    vector<double> reconstructed = computeIDFT(dft);

    //filtracja 1D
    vector<double> kernel = { 0.25, 0.5, 0.25 };
    vector<double> filtered = filter_1D(reconstructed, kernel);

    //Filtracja 2D
    // 2D Filtracja: u¿ycie y1 i y2 jako 2D sygna³u
    vector<vector<double>> signal2D = { y1, y2 };

    // 2D Gaussian blur kernel
    vector<vector<double>> kernel2D = {
        {0.0625, 0.125, 0.0625},
        {0.125,  0.25,  0.125},
        {0.0625, 0.125, 0.0625}
    };

    // Filtracja 2D
    vector<vector<double>> filtered2D = Filter2D(signal2D, kernel2D);

    // Wykresy
    subplot(3, 3, 1);
    plot(x, y1, "b-", x, y2, "g-");
    title("1i2");

    subplot(3, 3, 2);
    plot(x, y_combined);
    title("1+2");

    subplot(3, 3, 3);
    plotDFT(dft);

    subplot(3, 3, 4);
    plot(x, reconstructed);
    title("Rekonstrukcja z IDFT");

    subplot(3, 3, 5);
    plot(x, filtered);
    title("1D");

    subplot(3, 3, 6);
    imagesc(signal2D);
    title("bez 2D");

    subplot(3, 3, 7);
    imagesc(filtered2D);
    title("2D");
    show();

}