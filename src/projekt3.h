#pragma once

#ifndef PROJEKT3_H
#define PROJEKT3_H

#include <vector>
#include <string>

std::vector<double> operator*(const std::vector<double>& vec, double scalar);
std::vector<double> operator*(double scalar, const std::vector<double>& vec);

std::vector<double> generate_sine(double freq, double t_start, double t_end, size_t num_samples);
std::vector<double> generate_cosine(double freq, double t_start, double t_end, size_t num_samples);
std::vector<double> generate_square(double freq, double t_start, double t_end, size_t num_samples);
std::vector<double> generate_sawtooth(double freq, double t_start, double t_end, size_t num_samples);

struct DFTResult {
    std::vector<double> real;
    std::vector<double> imag;
};

DFTResult computeDFT(const std::vector<double>& signal);
std::vector<double> computeIDFT(const DFTResult& dft);

std::vector<double> filter_1D(const std::vector<double>& signal, const std::vector<double>& kernel);
std::vector<std::vector<double>> Filter2D(const std::vector<std::vector<double>>& signal, const std::vector<std::vector<double>>& filter);
std::vector<double> threshold_signal(const std::vector<double>& signal, double threshold);

void plotting(const std::string& type1, double freq1, double amp1,
    const std::string& type2, double freq2, double amp2,
    double t_start, double t_end, size_t num_samples,
    double threshold);

#endif // PROJEKT3_H