#include "HEAAN.h"

using namespace std;
using namespace NTL;

void Sqrt(Ciphertext &res, Ciphertext &x, int d, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    auto output = [&](Ciphertext& x, string name) {
        complex<double> sRes = scheme.decryptSingle(secretKey, x);
        cerr << "  " << name << " = " << sRes << ", logp = " << x.logp << ", logq = " << x.logq<< endl;
    };

    Ciphertext a0;
    a0.copy(x);
    Ciphertext b0;
    scheme.addConstAndEqual(x, -1., logp);  // x = x - 1
    b0.copy(x);

    Ciphertext a1;
    Ciphertext b;
    Ciphertext c;
    Ciphertext b1; 

    for (int n = 0; n < d; n++) {
        scheme.multByConst(b, b0, -0.5, logp);  // b = -b0 / 2
        scheme.reScaleByAndEqual(b, logp);
        scheme.addConstAndEqual(b, 1., logp);     // b = 1 + b      
        scheme.mult(a1, a0, b); // a1 = a0 * b
        scheme.modDownByAndEqual(a1, logp);
        scheme.reScaleByAndEqual(a1, logp);
        scheme.square(b, b0); // b = b0^2
        scheme.reScaleByAndEqual(b, logp);
        scheme.addConst(c, b0, -3, -1);  // c = b0 - 3
        scheme.divByPo2AndEqual(c, 2);  // c = c / 4
        scheme.mult(b1, b, c); // b1 = (b0^2) * ((b0 - 3) / 4)
        scheme.reScaleByAndEqual(b1, logp);

        a0.copy(a1);
        b0.copy(b1);
    }
    res.copy(a0);
}

void MinMax(Ciphertext &Min, Ciphertext &Max, Ciphertext &cipherA, Ciphertext &cipherB, int d, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    auto output = [&](Ciphertext& x, string name) {
        complex<double> sRes = scheme.decryptSingle(secretKey, x);
        cerr << "  " << name << " = " << sRes << ", logp = " << x.logp << ", logq = " << x.logq<< endl;
    };

    Ciphertext res;
    scheme.add(res, cipherA, cipherB);  // res = a + b

    Ciphertext x;
    scheme.multByConst(x, res, 0.5, logp);  // x = (a + b) / 2
    scheme.reScaleByAndEqual(x, logp);

    scheme.sub(res, cipherA, cipherB);   // res = a - b

    Ciphertext y;
    scheme.multByConst(y, res, 0.5, logp);  // y = (a - b) / 2
    scheme.reScaleByAndEqual(y, logp);

    Ciphertext z;
    scheme.mult(z, y, y);  // z = y^2
    scheme.reScaleByAndEqual(z, logp);

    Sqrt(res, z, d, logp, logq, scheme, secretKey);
    
    long bitsDown = x.logq - res.logq;
	scheme.modDownByAndEqual(x, bitsDown);

    scheme.sub(Min, x, res);  // min = x - res
    scheme.add(Max, x, res);  // max = x + res
}

void Inv(Ciphertext &res, Ciphertext &x, int d, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    auto output = [&](Ciphertext& x, string name) {
        complex<double> sRes = scheme.decryptSingle(secretKey, x);
        cerr << "  " << name << " = " << sRes << ", logp = " << x.logp << ", logq = " << x.logq<< endl;
    };

    Ciphertext a0;
    scheme.negateAndEqual(x);  // -x
    scheme.addConst(a0, x, 2., logp);  // x = 2 - x

    Ciphertext b0;
    scheme.addConst(b0, x, 1., logp);  // x = 1 - x

    Ciphertext a1;
    Ciphertext b;
    Ciphertext b1; 

    for (int n = 0; n < d; n++) {
        scheme.square(b1, b0); // b1 = b0^2
        scheme.reScaleByAndEqual(b1, logp);

        scheme.addConst(b, b1, 1., logp);  // b = b1 + 1
        scheme.mult(a1, a0, b); // a1 = a0 * b
        scheme.modDownByAndEqual(a1, logp);
        scheme.reScaleByAndEqual(a1, logp);

        a0.copy(a1);
        b0.copy(b1);
    }
    res.copy(a0);
}

void Comp(Ciphertext &result, Ciphertext &cipherA, Ciphertext &cipherB, int d, int d1, int t, long m, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    auto output = [&](Ciphertext& x, string name) {
        complex<double> sRes = scheme.decryptSingle(secretKey, x);
        cerr << "  " << name << " = " << sRes << ", logp = " << x.logp << ", logq = " << x.logq<< endl;
    };
    SchemeAlgo algo(scheme);

    Ciphertext res;
    scheme.add(res, cipherA, cipherB);  // res = a + b

    Ciphertext x;
    scheme.multByConst(x, res, 0.5, logp);  // x = (a + b) / 2
    scheme.reScaleByAndEqual(x, logp);

    Inv(res, x, d1, logp, logq, scheme, secretKey); // res = inv((a+b)/2)

    scheme.multByConst(x, cipherA, 0.5, logp);  // x = a / 2
    scheme.reScaleByAndEqual(x, logp);

    long bitsDown = x.logq - res.logq;
	scheme.modDownByAndEqual(x, bitsDown);

    Ciphertext a0;
    scheme.multAndEqual(x, res); // x = x * res
    scheme.reScaleByAndEqual(x, logp);
    a0.copy(x);  // a0 = x

    Ciphertext b0;
    scheme.negate(x, a0);  // x = -a0
    scheme.addConst(b0, x, 1., logp);  // b0 = 1 - a0

    Ciphertext a1;
    Ciphertext a;
    Ciphertext b;
    Ciphertext c;
    Ciphertext b1; 

    for (int n = 0; n < t; n++) {
        algo.power(a, a0, logp, m); // a = a0^m
        algo.power(b, b0, logp, m); // b = b0^m
        scheme.add(c, a, b);   // c = a + b
        Inv(res, c, d, logp, logq, scheme, secretKey); // res = inv(c)
        long bitsDown = a.logq - res.logq;
	    scheme.modDownByAndEqual(a, bitsDown);
        scheme.mult(a1, a, res); // a1 = a * res
        scheme.modDownByAndEqual(a1, logp);
        scheme.reScaleByAndEqual(a1, logp);

        scheme.negate(b1, a1);  // b1 = -a1
        scheme.addConstAndEqual(b1, 1., logp); // b1 = 1 + b1 = 1 - a1

        a0.copy(a1);
        b0.copy(b1);
    }
    result.copy(a0);
}

void ArrayMax(Ciphertext &result, Ciphertext *arr, int n, int d1, long logn, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    auto output = [&](Ciphertext& x, string name) {
        complex<double> sRes = scheme.decryptSingle(secretKey, x);
        cerr << "  " << name << " = " << sRes << ", logp = " << x.logp << ", logq = " << x.logq<< endl;
    };
    Ciphertext Min, Max;
    vector<vector<Ciphertext>> a(n + 1, vector<Ciphertext>(logn + 1, 0));
    for (int i = 1; i <= n; i++) {
        a[i][0].copy(arr[i - 1]);
    }
    int d = n;
    for (int j = 0; j < logn; j++) {
       cout << "here\n";
        if (d % 2 == 1) {
            a[d/2][j + 1].copy(a[d][j]);
        }
        d = n / 2;
        for (int i = 1; i <= d; i++) {
            Ciphertext ax, bx;
            ax.copy(a[2*i-1][j]);
            bx.copy(a[2*i][j]);
            MinMax(Min, Max, ax, bx, d1, logp, logq, scheme, secretKey);
            a[i][j + 1].copy(Max);
        }
    }
    result.copy(a[1][logn]);
}

Ciphertext* MaxIdx(Ciphertext *arr, int n, int d, int d1, int m, int t, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    auto output = [&](Ciphertext& x, string name) {
        complex<double> sRes = scheme.decryptSingle(secretKey, x);
        cerr << "  " << name << " = " << sRes << ", logp = " << x.logp << ", logq = " << x.logq<< endl;
    };
    SchemeAlgo algo(scheme);
    Ciphertext *b = new Ciphertext[n];
    Ciphertext inv, sum, x, y, sumb, sumbm;
    int first = 0;
    sum.copy(arr[0]);  // sum = arr[0]
    for (int i = 2; i <= n; i++) {
        scheme.addAndEqual(sum, arr[i - 1]);  // sum += arr[i-1]
    }

    double invN = 1.0/double(n);  // invN = 1/n
    scheme.multByConst(x, sum, invN, logp);  // x = sum/n
    scheme.reScaleByAndEqual(x, logp);

    Inv(inv, x, d1, logp, logq, scheme, secretKey);  // inv = inv(sum/n)

    for (int j = 1; j <= n - 1; j++) {
        scheme.multByConst(x, arr[j - 1], invN, logp); //x = arr[j-1] / n
        scheme.reScaleByAndEqual(x, logp);

        long bitsDown = x.logq - inv.logq;
	    scheme.modDownByAndEqual(x, bitsDown);
        scheme.mult(y, x, inv);  // y = arr[j-1] / n * inv
        scheme.reScaleByAndEqual(y, logp);

        b[j-1].copy(y);  // b[j - 1] = y
        if (first == 0) {
            sumb.copy(b[j-1]);
            first = 1;
            continue;
        }
        scheme.addAndEqual(sumb, b[j - 1]);  // sumb += b[j-1]
    }

    scheme.negate(x, sumb);  // x = -sumb
    scheme.addConst(b[n - 1], x, 1., logp);  // b[n-1] = 1 - sumb
    first = 0;
    for (int i = 1; i <= n; i++) {
        algo.power(y, b[i - 1], logp, m);  // y = b[i-1]^m
        if (first == 0) {
            sumbm.copy(y);
            first = 1;
            continue;
        }
        scheme.addAndEqual(sumbm, y);  // sumbm += b[i-1]^m
    }
    Inv(inv, sumbm, d, logp, logq, scheme, secretKey);  // inv = inv(sumbm)

    for (int i = 1; i <= t; i++) {
        for (int j = 0; j < n; j++) {
            algo.power(y, b[j], logp, m);  // y = b[j]^m
            scheme.mult(b[j], y, inv);  // b[j] = b[j]^m * inv
            scheme.reScaleByAndEqual(b[j], logp);
        }
        scheme.addConst(b[n - 1], x, 1., logp); 
    }
    return b;
}

void Threshold(Ciphertext &result, Ciphertext *arr, int n, int m, Ciphertext &b, int d, int d1, int t, long logn, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    auto output = [&](Ciphertext& x, string name) {
        complex<double> sRes = scheme.decryptSingle(secretKey, x);
        cerr << "  " << name << " = " << sRes << ", logp = " << x.logp << ", logq = " << x.logq<< endl;
    };
    Ciphertext* c = new Ciphertext[n];

    for (int i = 1; i <= n; i++) {
        Comp(c[i - 1], arr[i - 1], b, d, d1, t, m, logp, logq, scheme, secretKey);
    }
    Ciphertext sum;
    sum.copy(c[0]);  // sum = c[0]
    for (int i = 2; i <= n; i++) {
        scheme.addAndEqual(sum, c[i - 1]);  // sum += c[i-1]
    }
    result.copy(sum);
}

Ciphertext* TopKMax(Ciphertext *arr, int n, int m1, int d, int d1, int t, long logn, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    auto output = [&](Ciphertext& x, string name) {
        complex<double> sRes = scheme.decryptSingle(secretKey, x);
        cerr << "  " << name << " = " << sRes << ", logp = " << x.logp << ", logq = " << x.logq<< endl;
    };
    SchemeAlgo algo(scheme);
    Ciphertext *b = new Ciphertext[n];
    Ciphertext *c = new Ciphertext[n];
    Ciphertext *m = new Ciphertext[n];
    for (int i = 0; i < n; i++) {
        c[i].copy(arr[i]);
    }
    int first;
    Ciphertext prod, sum, negb, x;
    for (int j = 0; j < n; j++) {
        b = MaxIdx(c, n, d, d1, m1, t, logp, logq, scheme, secretKey);
        first = 0;
        for (int i = 0; i < n; i++) {
            scheme.mult(prod, b[i], c[i]);  // prod = b[i] * c[i]
            scheme.reScaleByAndEqual(prod, logp);
            if (first == 0) {
                sum.copy(prod);
                first = 1;
                continue;
            }
            scheme.addAndEqual(sum, prod);  // sum += b[i] * c[i]
        }
        m[j].copy(sum);
        for (int i = 0; i < n; i++) {
            scheme.negate(negb, b[i]);  // negb = -b[i]
            scheme.addConst(x, negb, 1., logp);  // x = 1 - b[i]
            scheme.mult(prod, c[i], x);  // prod = b[i] * c[i]
            scheme.reScaleByAndEqual(prod, logp);
            c[i].copy(prod);
        }
    }
    return m;
}

void TestSqrt(int d, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    TimeUtils timeutils;
    logq = 930;
    logp = 40;
    d = 11;

    double a = 0.5;

    Ciphertext cipherA;
    scheme.encryptSingle(cipherA, a, logp, logq);

    Ciphertext res;
    timeutils.start("Sqrt");
    Sqrt(res, cipherA, d, logp, logq, scheme, secretKey);
    timeutils.stop("Sqrt");

    complex<double> r = scheme.decryptSingle(secretKey, res);
    cout << "Sqrt result for " << a << " is: " << r << endl;
}

void TestMinMax(int d, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    TimeUtils timeutils;
    logq = 1170;
    logp = 30;
    d = 11;

    double a = 0.5;
    double b = 0.9;

    Ciphertext cipherA;
    scheme.encryptSingle(cipherA, a, logp, logq);
    
    Ciphertext cipherB;
    scheme.encryptSingle(cipherB, b, logp, logq);

    Ciphertext Min, Max;
    timeutils.start("MinMax");
    MinMax(Min, Max, cipherA, cipherB, d, logp, logq, scheme, secretKey);
    timeutils.stop("MinMax");

    complex<double> minRes = scheme.decryptSingle(secretKey, Min);
    complex<double> maxRes = scheme.decryptSingle(secretKey, Max);

    cout << "Min result is: " << minRes << endl;
    cout << "Max result is: " << maxRes << endl;
}

void TestInv(int d, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    TimeUtils timeutils;
    logq = 1600;
    logp = 40;
    d = 3;

    complex<double> a = 0.5;

    Ciphertext cipherA;
    scheme.encryptSingle(cipherA, a, logp, logq);

    Ciphertext res;
    timeutils.start("Inv");
    Inv(res, cipherA, d, logp, logq, scheme, secretKey);
    timeutils.stop("Inv");

    complex<double> r = scheme.decryptSingle(secretKey, res);
    cout << "Inv result for " << a << " is: " << r << endl;
}

void TestComp(int d, int d1, int t, int m, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    TimeUtils timeutils;
    logq = 2200;
    logp = 20;
    d = 5;
    d1 = 5;
    t = 6;
    m = 4;
    double a = 0.5;
    double b = 0.9;

    Ciphertext cipherA;
    scheme.encryptSingle(cipherA, a, logp, logq);
    
    Ciphertext cipherB;
    scheme.encryptSingle(cipherB, b, logp, logq);

    Ciphertext res;
    timeutils.start("Comp");
    Comp(res, cipherA, cipherB, d, d1, t, m, logp, logq, scheme, secretKey);
    timeutils.stop("Comp");

    complex<double> r = scheme.decryptSingle(secretKey, res);

    cout << "Comp result is: " << r << endl;
}

void TestArrayMax(int d, int n, long logn, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    TimeUtils timeutils;
    logq = 1600;
    logp = 40;
    logn = 1;
    n = 1 << logn;
    d = 8;
    vector<complex<double>> aux = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    complex<double>* vec = new complex<double>[n];
    for (int i = 0; i < n; i++) {
        vec[i] = aux[i % 9];
    }

    Ciphertext* arr = new Ciphertext[n];
    Ciphertext res;

    for (int i = 0; i < n; i++) {
        scheme.encryptSingle(arr[i], vec[i], logp, logq);
    }

    timeutils.start("ArrayMax");
    ArrayMax(res, arr, n, d, logn, logp, logq, scheme, secretKey);
    timeutils.stop("ArrayMax");

    complex<double> r = scheme.decryptSingle(secretKey, res);

    cout << "ArrayMax result is: " << r << endl;
}

void TestMaxIndx(int d, int d1, int n, int t, int m, long logn, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    TimeUtils timeutils;
    logq = 2050;
    logp = 20;
    logn = 1;
    n = 1 << logn;
    d = 11;
    d1 = 3;
    t = 3;
    m = 4;
    vector<complex<double>> aux = {0.5, 0.6, 0.7, 0.8, 0.9};
    complex<double>* vec = new complex<double>[n];
    for (int i = 0; i < n; i++) {
        vec[i] = aux[i % 5];
    }

    Ciphertext* arr = new Ciphertext[n];
    Ciphertext* res;

    for (int i = 0; i < n; i++) {
        scheme.encryptSingle(arr[i], vec[i], logp, logq);
    }
    
    timeutils.start("MaxIdx");
    res = MaxIdx(arr, n, d, d1, m, t, logp, logq, scheme, secretKey);
    timeutils.stop("MaxIdx");
    
    for (int i = 0; i < n; i++) {
        complex<double> r = scheme.decryptSingle(secretKey, res[i]);
        cout << "MaxIdx result is: " << r << endl;

    }
}

void TestThreshold(int d, int d1, int n, int t, int m, long logn, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    TimeUtils timeutils;
    logq = 2080;
    logp = 20;
    logn = 1;
    n = 1 << logn;
    d = 5;
    d1 = 5;
    t = 6;
    m = 4;
    complex<double> a = 0.5;
    Ciphertext cipherA;
    scheme.encryptSingle(cipherA, a, logp, logq);

    vector<complex<double>> aux = {0.5, 0.6, 0.7, 0.8, 0.9};
    complex<double>* vec = new complex<double>[n];
    for (int i = 0; i < n; i++) {
        vec[i] = aux[i % 5];
    }

    Ciphertext* arr = new Ciphertext[n];

    for (int i = 0; i < n; i++) {
        scheme.encryptSingle(arr[i], vec[i], logp, logq);
    }

    Ciphertext result;

    timeutils.start("Threshold");
    Threshold(result, arr, n, m, cipherA, d, d1, t, logn, logp, logq, scheme, secretKey);
    timeutils.stop("Threshold");
    complex<double> r = scheme.decryptSingle(secretKey, result);
    cout << "Threshold result is: " << r << endl;
}

void TestTopKMax(int d, int d1, int n, int t, int m, long logn, long logp, long logq, Scheme &scheme, SecretKey &secretKey) {
    TimeUtils timeutils;
    logq = 3091;
    logp = 40;
    logn = 1;
    n = 1 << logn;
    d = 4;
    d1 = 3;
    t = 3;
    m = 2;
    
    vector<complex<double>> aux = {0.5, 0.6, 0.7, 0.8, 0.9};
    complex<double>* vec = new complex<double>[n];
    for (int i = 0; i < n; i++) {
        vec[i] = aux[i % 5];
    }

    Ciphertext* arr = new Ciphertext[n];

    for (int i = 0; i < n; i++) {
        scheme.encryptSingle(arr[i], vec[i], logp, logq);
    }

    Ciphertext* res;

    timeutils.start("TopKMax");
    res = TopKMax(arr, n, m, d, d1, t, logn, logp, logq, scheme, secretKey);
    timeutils.stop("TopKMax");
 
    for (int i = 0; i < n; i++) {
        complex<double> r = scheme.decryptSingle(secretKey, res[i]);
        cout << "TopKMax result is: " << r << endl;
    }
}

int main() {

    long logq = 1600; ///< Ciphertext modulus (this value should be <= logQ in "scr/Params.h")
    long logp = 40; ///< Scaling Factor (larger logp will give you more accurate value)
    long logn = 1; ///< number of slot is 1024 (this value should be < logN in "src/Params.h")
    long n = 1 << logn;
    long slots = n;
    long numThread = 8;
    int d = 3;
    int d1 = 3;
    int t = 3;
    long m = 2;

    srand(time(NULL));
    SetNumThreads(numThread);
    TimeUtils timeutils;
    Ring ring;
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);
    //scheme.addMultKey(secretKey);
    //scheme.addLeftRotKeys(secretKey); ///< When you need left rotation for the vectorized message
    //scheme.addRightRotKeys(secretKey); ///< When you need right rotation for the vectorized message

    TestSqrt(d, logp, logq, scheme, secretKey);
    TestMinMax(d, logp, logq, scheme, secretKey);
    TestInv(d, logp, logq, scheme, secretKey);
    // TestComp(d, d1, t, m, logp, logq, scheme, secretKey);
    // TestArrayMax(d, n, logn, logp, logq, scheme, secretKey);
    // TestMaxIndx(d, d1, n, t, m, logn, logp, logq, scheme, secretKey);
    // TestThreshold(d, d1, n, t, m, logn, logp, logq, scheme, secretKey);
    // TestTopKMax(d, d1, n, t, m, logn, logp, logq, scheme, secretKey);

    return 0;
}