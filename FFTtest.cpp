#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <algorithm>
#include <iomanip>

using namespace std;

class MixedRadixFFT {
private:
    using Complex = complex<double>;
    using ComplexVector = vector<Complex>;
    
    const double PI = acos(-1.0);
        
    bool isSmoothNumber(int n) {
        if (n <= 0) return false;
        while (n % 2 == 0) n /= 2;
        while (n % 3 == 0) n /= 3;
        while (n % 5 == 0) n /= 5;
        return n == 1;
    }
    
    vector<int> getPrimeFactors(int n) {
        vector<int> factors;
        for (int factor : {2, 3, 5}) {
            while (n % factor == 0) {
                factors.push_back(factor);
                n /= factor;
            }
        }
        return factors;
    }
    
public:
    
    ComplexVector forward(const ComplexVector& input) {
        int n = input.size();
        if (!isSmoothNumber(n)) {
            throw invalid_argument("Длина преобразования должна быть кратной только 2, 3 и 5");
        }
        
        ComplexVector result = input;
        vector<int> factors = getPrimeFactors(n);
        
        int offset = 1;
        for (int stage = 0; stage < (int)factors.size(); stage++) {
            int radix = factors[stage];
            int length = radix * offset;
            
            for (int start = 0; start < n; start += length) {
                for (int i = 0; i < offset; i++) {
                    
                    vector<Complex> temp(radix);
                    for (int r = 0; r < radix; r++) {
                        temp[r] = result[start + i + r * offset];
                    }
                    
                    
                    for (int r = 0; r < radix; r++) {
                        Complex sum = 0;
                        for (int k = 0; k < radix; k++) {
                            double angle = -2.0 * PI * r * k / radix;
                            Complex twiddle = Complex(cos(angle), sin(angle));
                            sum += twiddle * temp[k];
                        }
                        result[start + i + r * offset] = sum;
                    }
                }
            }
            offset = length;
        }
        
        
        if (all_of(factors.begin(), factors.end(), [](int f) { return f == 2; })) {
            int log2n = log2(n);
            ComplexVector output(n);
            for (int i = 0; i < n; i++) {
                int rev = 0;
                for (int j = 0; j < log2n; j++) {
                    if (i & (1 << j)) {
                        rev |= 1 << (log2n - 1 - j);
                    }
                }
                output[rev] = result[i];
            }
            return output;
        }
        
        return result;
    }
    
    
    ComplexVector inverse(const ComplexVector& input) {
        int n = input.size();
        ComplexVector result(n);
                
        for (int i = 0; i < n; i++) {
            result[i] = conj(input[i]);
        }
                
        result = forward(result);
                
        for (int i = 0; i < n; i++) {
            result[i] = conj(result[i]) / Complex(n, 0);
        }
        
        return result;
    }
    
    double computeRelativeError(const ComplexVector& original, const ComplexVector& reconstructed) {
        if (original.size() != reconstructed.size()) {
            throw invalid_argument("Векторы должны быть одинаковой длины");
        }
        
        double normOriginal = 0.0;
        double error = 0.0;
        
        for (size_t i = 0; i < original.size(); i++) {
            normOriginal += norm(original[i]);
            error += norm(original[i] - reconstructed[i]);
        }
        
        if (normOriginal < 1e-30) return error;
        return sqrt(error / normOriginal);
    }
        
    double computeMaxError(const ComplexVector& original, const ComplexVector& reconstructed) {
        if (original.size() != reconstructed.size()) {
            throw invalid_argument("Векторы должны быть одинаковой длины");
        }
        
        double maxError = 0.0;
        for (size_t i = 0; i < original.size(); i++) {
            maxError = max(maxError, abs(original[i] - reconstructed[i]));
        }
        return maxError;
    }
};

vector<complex<double>> generateRandomComplexVector(int n, double minVal = -10.0, double maxVal = 10.0) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(minVal, maxVal);
    
    vector<complex<double>> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = complex<double>(dist(gen), dist(gen));
    }
    return result;
}

void printVector(const vector<complex<double>>& vec, const string& name, int maxItems = 10) {
    cout << name << " [" << vec.size() << "]:\n";
    for (int i = 0; i < min(maxItems, (int)vec.size()); i++) {
        cout << fixed << setprecision(6);
        cout << "  [" << setw(2) << i << "] = " 
             << setw(10) << vec[i].real() << " + " 
             << setw(10) << vec[i].imag() << "i\n";
    }
    if (vec.size() > maxItems) {
        cout << "  ... (показано " << maxItems << " из " << vec.size() << ")\n";
    }
}

int main() {
    
    const int N = 30;  // Длина преобразования (должна быть кратной 2, 3, 5)
    const bool PRINT_DETAILS = true;  
    
    
    
    cout << "   Mixed-Radix FFT (2, 3, 5)\n";
    
    
    cout << "Длина преобразования: N = " << N << "\n";
        
    MixedRadixFFT fft;
    if (N <= 0 || !(N == 1 || (N % 2 == 0 || N % 3 == 0 || N % 5 == 0))) {
        cout << "Предупреждение: N = " << N << " может не разлагаться на множители 2, 3, 5\n";
    }
    
    try {
        
        vector<complex<double>> input = generateRandomComplexVector(N, -5.0, 5.0);
        
        if (PRINT_DETAILS) {
            printVector(input, "Исходный сигнал");
            cout << endl;
        }
                
        vector<complex<double>> freqDomain = fft.forward(input);
        
        if (PRINT_DETAILS) {
            printVector(freqDomain, "Частотная область");
            cout << endl;
        }
                
        vector<complex<double>> output = fft.inverse(freqDomain);
        
        if (PRINT_DETAILS) {
            printVector(output, "Восстановленный сигнал");
            cout << endl;
        }
                
        double relativeError = fft.computeRelativeError(input, output);
        double maxError = fft.computeMaxError(input, output);
        
        
        cout << "   Результаты\n";
        
        cout << "Относительная ошибка: " << scientific << setprecision(6) << relativeError << "\n";
        cout << "Максимальная ошибка:  " << scientific << maxError << "\n";
        
        if (relativeError < 1e-14) {
            
            cout << "Ошибка находится в пределах машинной точности.\n";
        } else {
            cout << "\nОшибка превышает ожидаемую.\n";
        }
        
        
        
    } catch (const exception& e) {
        cerr << "\nОШИБКА: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}