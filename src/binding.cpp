#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "projekt3.h"

namespace py = pybind11;

PYBIND11_MODULE(signals, m) {
    m.doc() = "Modul do generowania i analizy sygnałow";

    m.def("sin", &generate_sine,
        py::arg("freq"),
        py::arg("t_start"),
        py::arg("t_end"),
        py::arg("num_samples"),
        "Generuje sygnal sinusoidalny");

    m.def("cos", &generate_cosine,
        py::arg("freq"),
        py::arg("t_start"),
        py::arg("t_end"),
        py::arg("num_samples"),
        "Generuje sygnal cosinusoidalny");

    m.def("square", &generate_square,
        py::arg("freq"),
        py::arg("t_start"),
        py::arg("t_end"),
        py::arg("num_samples"),
        "Generuje sygnal prostokatny");

    m.def("sawtooth", &generate_sawtooth,
        py::arg("freq"),
        py::arg("t_start"),
        py::arg("t_end"),
        py::arg("num_samples"),
        "Generuje sygnal piloksztaltny");

    py::class_<DFTResult>(m, "DFTResult")
        .def_readwrite("real", &DFTResult::real)
        .def_readwrite("imag", &DFTResult::imag);

    m.def("dft", &computeDFT,
        py::arg("signal"),
        "Liczenie DFT z podanego sygnalu");

    m.def("idft", &computeIDFT,
        py::arg("dft_result"),
        "Odwrotna transformata Fouriera");

    m.def("filter_1D", &filter_1D,
        py::arg("signal"),
        py::arg("kernel"),
        "Filtracja sygnalu 1D za pomoca jądra filtrującego");

    m.def("threshold", &threshold_signal,
        py::arg("signal"),
        py::arg("threshold"),
        "Zwraca sygnal progowany: wartosci > prog to 1.0, reszta 0.0");

    m.def("filter_2D", &Filter2D,
        py::arg("signal"),
        py::arg("filter"),
        "Filtracja sygnalu 2D za pomoca jądra filtrującego");


    m.def("run", [](const std::string& type1, double freq1, double amp1,
        const std::string& type2, double freq2, double amp2,
        double t_start, double t_end, size_t num_samples,
        double threshold) {
            return plotting(type1, freq1, amp1, type2, freq2, amp2, t_start, t_end, num_samples, threshold);
        });


}