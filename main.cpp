#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <bitset>
#include <algorithm>

using namespace std;
int n, k;
vector<vector<int>> G_etalon;
vector<vector<int>> G;
string path = R"(C:\Users\S3\CLionProjects\order-statistic\source\RM_16_11_4.gen)";
//string path = R"(C:\Users\S3\CLionProjects\order-statistic\source\bch_32_21_6.gen)";

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

int rnd(int l, int r) {
    uniform_int_distribution<int> range(l, r);
    return range(rng);
}

std::random_device rd;
std::mt19937 e2(rd());

void getG() {
    ifstream in;
    in.open(path);
    in >> n >> k;
    G.resize(n);
    for (int i = 0; i < n; ++i) {
        G[i].resize(k);
    }

    for (int i = 0; i < k; ++i) {
        string s;
        in >> s;
        for (int j = 0; j < n; ++j) {
            G[j][i] = (s[j] - '0');
        }
    }
    G_etalon = G;
    in.close();
}

int mul(vector<int> &to, vector<int> &y) {
    if (to.size() != y.size()) {
        cout << "Incorrect size in mul";
        terminate();
    }
    int ans = 0;
    for (int i = 0; i < to.size(); ++i) {
        ans = (ans + to[i] * y[i]) % 2;
    }
    return ans;
}

vector<int> matrixMul(vector<int> &v) {
    if (k != v.size()) {
        cout << "Incorrect size in matrixMul";
        terminate();
    }

    vector<int> ans(n);
    for (int i = 0; i < n; ++i) {
        ans[i] = mul(v, G[i]);
    }
    return ans;
}

vector<int> get_random_codeword() {
    int data = rnd(0, (1 << k) - 1);
    bitset<32> bitset(data);
    vector<int> data_vector(k);
    for (int i = 0; i < k; ++i) {
        data_vector[i] = bitset[i];
    }
    return matrixMul(data_vector);
}

vector<double> get_corrupt(vector<int> &random_codeword, double arg) {
    std::normal_distribution<double> dist(0, sqrt(static_cast<double>(n) / (2 * arg * k)));
    vector<double> ans(n, 0);
    for (int i = 0; i < n; ++i) {
        ans[i] = ((random_codeword[i] * 2) - 1) + dist(e2);
    }
    return ans;
}

vector<int> get_signs(vector<double> &corrupted) {
    vector<int> ans(n);
    for (int i = 0; i < n; ++i) {
        ans[i] = (corrupted[i] > 0) ? 1 : -1;
    }
    return ans;
}

vector<double> get_reliability(vector<double> &corrupted) {
    vector<double> ans(n);
    for (int i = 0; i < n; ++i) {
        ans[i] = abs(corrupted[i]);
    }
    return ans;
}


vector<int> get_permutation(vector<double> &reliability) {
    std::vector<pair<double, int>> elements;
    for (int i = 0; i < reliability.size(); ++i) {
        elements.emplace_back(reliability[i], i);
    }
    auto comparator =
            [](const pair<double, int> &a, const pair<double, int> &b) { return a.first > b.first; };
    std::sort(elements.begin(), elements.end(), comparator);
    std::vector<int> result;
    result.reserve(reliability.size());
    for (const auto &element: elements) {
        result.push_back(element.second);
    }
    return result;
}

void permuteG(vector<int> &permutation) {
    for (int i = 0; i < n; ++i) {
        G[i] = G_etalon[permutation[i]];
    }
}

void sumLines(int to, int from) {
    for (int i = 0; i < n; ++i) {
        G[i][to] += G[i][from];
        G[i][to] %= 2;
    }
}

void gauss(vector<int> &permutatin) {
    for (int i = 0; i < k; ++i) {
        if (G[i][i] == 0) {
            int j;
            for (j = i + 1; j < k; ++j) {
                if (G[i][j] != 0) {
                    break;
                }
            }
            if (j != k) {
                sumLines(i, j);
            } else {
                int l;
                for (l = i + 1; l < n; ++l) {
                    if (G[l][i] == 1) {
                        break;
                    }
                }
                swap(G[i], G[l]);
                swap(permutatin[i], permutatin[l]);
            }
        }

        for (int j = i + 1; j < k; ++j) {
            if (G[i][j] != 0) {
                sumLines(j, i);
            }
        }
    }
    for (int i = k - 1; i > 0; --i) {
        for (int j = i - 1; j >= 0; --j) {
            if (G[i][j] != 0) {
                sumLines(j, i);
            }
        }
    }

}


vector<int> order_statistic(int t, vector<int>& signs, vector<int>& permutation) {
    vector<int> y(k);
    for (int i = 0; i < k; ++i) {
        y[i] = (signs[permutation[i]] + 1) / 2;
    }
    vector<int> ans;
    long double mu_min = numeric_limits<long double>::max();
    for (int i = 0; i < (1 << k); ++i) {
        bitset<32> bitset(i);
        if (bitset.count() <= t){
            vector<int> withError(k);
            for (int j = 0; j < k; ++j) {
                withError[j] = (bitset[j] + y[j]) % 2;
            }

            vector<int> result = matrixMul(withError);
            for (int j = 0; j < n; ++j) {
                result[j] = 2 * result[j] - 1;
            }

            long double mu = 0;
            for (int j = 0; j < n; ++j) {
                mu += (result[j] - signs[permutation[j]]) * (result[j] - signs[permutation[j]]);
            }
            mu = sqrt(mu);
            if (mu < mu_min) {
                mu_min = mu;
                ans = result;

            }
        }
    }
    return ans;
}

int main() {
    getG();
    for (double Eb = 3.0; Eb < 6.0; Eb += 0.1) {
        double correct = 0;
        double all = 3000;
        for(int tries = 0; tries < all; tries++){
            G = G_etalon;
            vector<int> random_codeword = get_random_codeword();
            vector<double> corrupt = get_corrupt(random_codeword, Eb);
            vector<int> signs = get_signs(corrupt);
            vector<double> reliability = get_reliability(corrupt);
            vector<int> permutation = get_permutation(reliability);
            permuteG(permutation);
            gauss(permutation);
            vector<int> result = order_statistic(2, signs, permutation);
            for (int i = 0; i < n; ++i) {
                result[i] = (result[i] + 1) / 2;
            }
            bool flag = true;
            int count_true = 0;
            int count_res = 0;
            for (int i = 0; i < n; ++i) {
                if (result[i] != random_codeword[permutation[i]]) {
                    flag = false;
                }
                if (result[i] == 1) {
                    count_res++;
                }
                if (random_codeword[i] == 1) {
                    count_true++;
                }
            }

            if (flag) {
                correct++;
            }
        }
        cout << fixed << setprecision(5) << Eb << ";" << correct / all << '\n';
    }
}
