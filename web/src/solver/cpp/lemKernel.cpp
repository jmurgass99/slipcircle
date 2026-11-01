// lemKernel.cpp — C++ WebAssembly kernel for Bishop Simplified m_alpha and
// seismic driving moment evaluation.

extern "C" {

double bishopMAlpha(double cos_a, double sin_a, double tan_phi, double fs) {
    return cos_a + (sin_a * tan_phi) / fs;
}

double sliceDrivingTerm(double w, double sin_a, double kh, double arm) {
    return w * sin_a + w * kh * arm;
}

}
