#include <iostream>
#include "BigInteger.h"
#include "Algorithms.h"
#include <map>
using  namespace std;
using namespace LongArithmetic;

enum Algorithm{
    PollardRho,
    BabyStepGiantStep,
    EulerFunc,
    MobiusFunc,
    Legendre,
    Jacobi,
    Cipolla,
    MillerRabin,
    ElGamal
};

static map<int, Algorithm> int_to_algo =
        {
                {1, PollardRho},
                {2, BabyStepGiantStep},
                {3, EulerFunc},
                {4, MobiusFunc},
                {5, Legendre},
                {6, Jacobi},
                {7, Cipolla},
                {8, MillerRabin},
                {9,ElGamal}
        };

int main() {
    Algorithms solve = Algorithms();
    EllipticCurve curve = EllipticCurve();

    int algo;
    BigInteger n; // value for factorization, prime, euler, mobius ; message for el gamal
    BigInteger a, b, m; // for discrete log
    BigInteger A, p; // for legendre and jacobi; for discrete root
    pair<BigInteger, BigInteger> res; // for discrete root
    while(true) {
        cout << "\nChoose action:\n"
                "1 - Pollard Rho Factorization\n"
                "2 - Discrete Logarithm (Baby-step-Giant-step)\n"
                "3 - Euler Function\n"
                "4 - Mobius Function\n"
                "5 - Legendre Symbol\n"
                "6 - Jacobi Symbol\n"
                "7 - Discrete Square Root\n"
                "8 - Is Prime (Miller-Rabin)\n"
                "9 - El-Gamal encryption over elliptic curve\n";
        cin >> algo;
        if (int_to_algo[algo] == PollardRho || int_to_algo[algo] == EulerFunc || int_to_algo[algo] == MobiusFunc || int_to_algo[algo] == MillerRabin) {
            cout << "Enter a number:\n";
            cin >> n;
        }

        switch(int_to_algo[algo]) {
            case PollardRho:
                if (solve.is_Prime(n)) {
                    cout << n << " is prime.";
                    continue;
                }
                else {
                    //auto factors = solve.pollardRhoFactorization(n);
                    auto factors = solve.pollardRho(n);
                    cout << "Factorization for " << n << " = " << factors;
//                    for (auto f : factors) {
//                        cout << f << " ";
//                    }
                    //cout << "\n";
                }
                break;
            case BabyStepGiantStep:
                cout << "The congruence equation is a^x = b mod m, where a - log base, b - , m - modulo, x to be found.\n"
                        "Enter a, b and m:\n";
                cin >> a >> b >> m;
//                if (gcd(a, m) > 1) {
//                    cout << "Numbers a and m are not co-prime.\n";
//                    continue;
//                }
//                else {}
                cout << "x = " << solve.discreteLogarithm(a, b, m);
                break;
            case EulerFunc:
                cout << "phi( " << n << " ) = " << solve.eulerFunction(n);
                break;
            case MobiusFunc:
                cout << "mu( " << n << " ) = " << solve.mobiusFunction(n);
                break;
            case Legendre:
                cout << "Legendre symbol a / p, where a - big integer, p - prime number.\n"
                        "Enter a and p:\n";
                cin >> A >> p;
                if (!solve.is_Prime(p)) {
                    cout << "p must be prime.";
                    continue;
                }
                else {
                    cout << "L(a, p) = " << solve.legendreSymbol(A, p);
                }
                break;
            case Jacobi:
                cout << "Enter a, n for Jacobi symbol (a / n):\n";
                cin >> A >> p;
                cout << "J(a, p) = " << solve.jacobiSymbol(A, p);
                break;
            case Cipolla:
                cout << "The congruence equation is x^2 = n mod p, x to be found.\n"
                        "Enter n and p:\n";
                cin >> A >> p;
                cout << "x^2 = " << A << " mod " << p;
//                if (!solve.is_Prime(p)) {
//                    cout << "p must be prime.";
//                    continue;
//                }
//                else {
                    res = solve.discreteSquareRoot(A, p);
                    cout << "\nx1 = " << res.first << "\nx2 = " << res.second;
                //}
                break;
            case MillerRabin:
                if(solve.is_Prime(n)){
                    cout << n << " is prime.";
                }
                else
                    cout << n << " is not a prime number.";
                break;
            case ElGamal:
                cout << "Enter the message (a number only)\n";
                cin >> n;
                cout << "Message is: " << n;
                cout << "\n Decrypted message: " << solve.decryption(curve, n);
                break;
            default:
                cout << "There is no such action";
                break;
        }
    }

    return 0;
}