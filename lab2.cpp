#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;

//======================================================================================================================
//                                                   CLASS Poly
//======================================================================================================================

class Poly {
public:
    int p;                  // модуль (простое число)
    std::vector<int> a;     // коэффициенты: a[0] — свободный, a[k] — x^k

public:
    Poly(int mod = 2)
        : p(mod){
        a.clear();
        a.push_back(0);
    }

    Poly(const std::vector<int>& v, int mod)
        : p(mod){
        a = v;
        Normalize();
    }

    Poly(const Poly& x)
        : p(x.p), a(x.a){}

    //==================================================== Normalize ====================================================

    void Normalize(){
        while (a.size() > 1 && a.back() % p == 0){
            a.pop_back();
        }

        for (int &x : a){
            x = (x % p + p) % p;
        }
    }

    int Degree() const {
        return (int)a.size() - 1;
    }

    //============================================================ Операторы ============================================================

    Poly operator + (const Poly& t) const {
        Poly r(*this);
        if (t.a.size() > r.a.size()){
            r.a.resize(t.a.size(), 0);
        }
        for (size_t i = 0; i < t.a.size(); i++){
            r.a[i] = (r.a[i] + t.a[i]) % p;
        }
        r.Normalize();
        return r;
    }

    Poly operator - (const Poly& t) const {
        Poly r(*this);
        if (t.a.size() > r.a.size()){
            r.a.resize(t.a.size(), 0);
        }
        for (size_t i = 0; i < t.a.size(); i++){
            r.a[i] = (r.a[i] - t.a[i] + p) % p;
        }
        r.Normalize();
        return r;
    }

    Poly operator * (const Poly& t) const {
        Poly r(p);
        r.a.assign(a.size() + t.a.size() - 1, 0);

        for (size_t i = 0; i < a.size(); i++){
            for (size_t j = 0; j < t.a.size(); j++){
                r.a[i + j] = (r.a[i + j] + a[i] * 1LL * t.a[j]) % p;
            }
        }

        r.Normalize();
        return r;
    }

    //====================================================== Деление с остатком =======================================================

    std::pair<Poly, Poly> DivMod(const Poly& d) const {
        if (d.Degree() == 0 && d.a[0] == 0){
            throw std::runtime_error("division by zero polynomial");
        }

        Poly R(*this);
        Poly Q(p);

        Q.a.assign(std::max(0, Degree() - d.Degree() + 1), 0);

        int invLead = InvMod(d.a.back(), p);

        while (R.Degree() >= d.Degree() && !(R.a.size() == 1 && R.a[0] == 0)){
            int coef = (R.a.back() * 1LL * invLead) % p;
            int shift = R.Degree() - d.Degree();

            Q.a[shift] = coef;

            std::vector<int> sub(shift + d.a.size());
            for (size_t i = 0; i < d.a.size(); i++){
                sub[shift + i] = (coef * 1LL * d.a[i]) % p;
            }

            R = R - Poly(sub, p);
        }

        Q.Normalize();
        R.Normalize();
        return {Q, R};
    }

    //=========================================================== GCD ===========================================================

    static Poly GCD(Poly A, Poly B){
        while (!(B.a.size() == 1 && B.a[0] == 0)){
            auto R = A.DivMod(B).second;
            A = B;
            B = R;
        }
        A.Normalize();
        return A;
    }

    //========================================= x^k mod f(x) =========================================

    static Poly PowX(long long k, const Poly& mod){
        Poly r({0,1}, mod.p); // r = x
        Poly res({1}, mod.p); // res = 1

        while (k > 0){
            if (k & 1){
                res = (res * r).DivMod(mod).second;
            }
            r = (r * r).DivMod(mod).second;
            k >>= 1;
        }
        return res;
    }

    //================================================= Проверка неприводимости =================================================

    static bool IsIrreducible(const Poly& f){
        int p = f.p;
        int n = f.Degree();

        if (n <= 0){
            return false;
        }

        Poly test = PowX(PowInt(p, n), f) - Poly({0,1}, p);
        if (!(test.a.size() == 1 && test.a[0] == 0)){
            return false;
        }

        std::vector<int> divs = PrimeDivisors(n);

        for (int q : divs){
            long long pw = PowInt(p, n / q);
            Poly h = PowX(pw, f) - Poly({0,1}, p);

            if (GCD(f, h).Degree() > 0){
                return false;
            }
        }

        return true;
    }

    //=========================================================== Вывод ===========================================================

    friend std::ostream& operator<<(std::ostream& out, const Poly& x){
        for (int i = x.Degree(); i >= 0; i--){
            out << x.a[i];
            if (i > 0){
                out << "x^" << i << " + ";
            }
        }
        return out;
    }

private:

    static int InvMod(int a, int p){
        int r0 = a, r1 = p;
        int s0 = 1, s1 = 0;

        while (r1){
            int q = r0 / r1;
            r0 -= q * r1; std::swap(r0, r1);
            s0 -= q * s1; std::swap(s0, s1);
        }

        if (r0 != 1){
            throw std::runtime_error("element is not invertible");
        }

        return (s0 % p + p) % p;
    }

    static long long PowInt(long long a, long long b){
        long long r = 1;
        while (b--){
            r *= a;
        }
        return r;
    }

    static std::vector<int> PrimeDivisors(int n){
        std::vector<int> d;
        for (int i = 2; i * i <= n; i++){
            if (n % i == 0){
                d.push_back(i);
                while (n % i == 0){
                    n /= i;
                }
            }
        }
        if (n > 1){
            d.push_back(n);
        }
        return d;
    }
};

//======================================================================================================================
//                                                    MAIN
//======================================================================================================================

Poly ReadPolyManual(int p){
    cout << "Enter the degree n: ";
    int n;
    cin >> n;

    vector<int> c(n + 1);
    for (int i = 0; i <= n; i++){
        cout << "The coefficient at x^" << i << ": ";
        cin >> c[i];
    }

    Poly f(c, p);
    cout << "The polynomial is introduced: " << f << endl;
    return f;
}

int main(){
    int p;

    cout << "Enter the modulus p (a prime number): ";
    cin >> p;

    cout << "\nEnter the polynomial f(x):" << endl;
    Poly f = ReadPolyManual(p);

    bool irr = Poly::IsIrreducible(f);

    cout << "Status: " << (irr ? "nepevodim" : "pevodim") << endl;

    return 0;
}

