#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include "IO.hpp"
#include "node.hpp"
#include "Elemstiffness.hpp"
#include "B_matrix.hpp"
#include "ja.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;
int main(int argc, char **argv)
{

    if (argc == 1)
    {
        throw invalid_argument("please  input file!");
    }
    auto d = new node(argv[1]);

    //初始化总体刚度矩阵
    vector<vector<double>>
        K(2 * d->nodes.size(), vector<double>(2 * d->nodes.size(), 0));
    //初始化单元刚度矩阵
    vector<vector<double>> Ke(8, vector<double>(8, 0));

    for (size_t i = 0; i < d->elements.size(); ++i)
    {
        Ke = ElementStiffness(166667, 0.3, 1, d->nodes[d->elements[i][0] - 1][0], d->nodes[d->elements[i][0] - 1][1],
                              d->nodes[d->elements[i][1] - 1][0], d->nodes[d->elements[i][1] - 1][1], d->nodes[d->elements[i][2] - 1][0],
                              d->nodes[d->elements[i][2] - 1][1], d->nodes[d->elements[i][3] - 1][0], d->nodes[d->elements[i][3] - 1][1]);

        for (int j = 0; j < 4; j++)
        {
            for (int p = 0; p < 4; p++)
            {
                K[2 * d->elements[i][j] - 2][2 * d->elements[i][p] - 2] = Ke[2 * j][2 * p] + K[2 * d->elements[i][j] - 2][2 * d->elements[i][p] - 2];
                K[2 * d->elements[i][j] - 2][2 * d->elements[i][p] - 2 + 1] = Ke[2 * j][2 * p + 1] + K[2 * d->elements[i][j] - 2][2 * d->elements[i][p] - 2 + 1];
                K[2 * d->elements[i][j] - 2 + 1][2 * d->elements[i][p] - 2] = Ke[2 * j + 1][2 * p] + K[2 * d->elements[i][j] - 2 + 1][2 * d->elements[i][p] - 2];
                K[2 * d->elements[i][j] - 2 + 1][2 * d->elements[i][p] - 2 + 1] = Ke[2 * j + 1][2 * p + 1] + K[2 * d->elements[i][j] - 2 + 1][2 * d->elements[i][p] - 2 + 1];
            }
        }
    }

    cout << "Stiffness calculate complete" << endl;

    for (size_t i = 0; i < d->bounds.size(); i++)
    {
        if (d->bounds[i][0] == 0)
            ; //不处理
        else
        {
            K[2 * i][2 * i] = 1;
            d->loads[i][0] = 0;
            // 一整行置0
            for (size_t j = 0; j < 2 * i; j++)
                K[2 * i][j] = 0;
            for (size_t j = 2 * i + 1; j < 2 * d->bounds.size(); j++)
                K[2 * i][j] = 0;

            for (size_t j = 0; j < 2 * i; j++)
                K[j][2 * i] = 0;
            for (size_t j = 2 * i + 1; j < 2 * d->bounds.size(); j++)
                K[j][2 * i] = 0;
        }
        if (d->bounds[i][1] == 0)
            ; //不处理
        else
        {
            K[2 * i + 1][2 * i + 1] = 1;
            d->loads[i][1] = 0;
            for (size_t j = 0; j < 2 * i + 1; j++)
                K[2 * i + 1][j] = 0;
            for (size_t j = 2 * i + 2; j < 2 * d->bounds.size(); j++)
                K[2 * i + 1][j] = 0;

            for (size_t j = 0; j < 2 * i + 1; j++)
                K[j][2 * i + 1] = 0;
            for (size_t j = 2 * i + 2; j < 2 * d->bounds.size(); j++)
                K[j][2 * i + 1] = 0;
        }
    }

    vector<double> F(2 * d->nodes.size(), 0);

    for (size_t i = 0; i < d->nodes.size(); i++) //把两列化为一列
    {
        F[2 * i] = d->loads[i][0];
        F[2 * i + 1] = d->loads[i][1];
    }
    cout << "Use eigen to solve equation" << endl;

    //用eigen求解方程
    MatrixXd Kmatrix(2 * d->nodes.size(), 2 * d->nodes.size());
    VectorXd Fmatrix(2 * d->nodes.size());

    for (size_t j = 0; j < 2 * d->nodes.size(); j++) //共2 * sizeof(nodes)行
    {
        for (size_t i = 0; i < 2 * d->nodes.size(); i++) //共2 * sizeof(nodes)列 组成一行
        {
            Kmatrix(j, i) = K[j][i];
        }
    }
    for (size_t i = 0; i < 2 * d->nodes.size(); i++) //共2 * sizeof(nodes)列 组成一行
    {
        Fmatrix(i) = F[i];
    }

    // QR 分解求出位移向量x
    VectorXd x = Kmatrix.colPivHouseholderQr().solve(Fmatrix);
    //把位移向量化为n乘2的矩阵，方便输出
    vector<vector<double>> X(d->nodes.size(), vector<double>(2, 0.0));
    for (size_t i = 0; i < 2 * d->nodes.size(); i = i + 2) //把一列化为两列
    {
        X[i / 2][0] = x(i);
        X[i / 2][1] = x(i + 1);
    }

    cout << "eigen over!" << endl;

    //应力
    vector<double> stress_x(d->nodes.size(), 0);
    vector<double> stress_y(d->nodes.size(), 0);
    vector<double> stress_xy(d->nodes.size(), 0);

    double E = 166667, u = 0.3;
    double D[3][3] = {E / (1 - u * u), E * u / (1 - u * u), 0, E * u / (1 - u * u), E / (1 - u * u), 0, 0, 0, (1 - u) * E / (2 * (1 - u * u))};
    double f[2] = {-0.577350, 0.577350}; //假设需要两个高斯积分点，查询高斯积分表得到两个高斯点的函数值
    for (size_t i = 0; i < d->elements.size(); i++)
    {
        vector<vector<double>> bmatrix = BMatrix(d->nodes[d->elements[i][0] - 1][0], d->nodes[d->elements[i][0] - 1][1],
                                                 d->nodes[d->elements[i][1] - 1][0], d->nodes[d->elements[i][1] - 1][1],
                                                 d->nodes[d->elements[i][2] - 1][0], d->nodes[d->elements[i][2] - 1][1],
                                                 d->nodes[d->elements[i][3] - 1][0], d->nodes[d->elements[i][3] - 1][1], f[0], f[1]);
        for (int j = 0; j < 4; j++)
        {
            double BX[3];
            BX[0] = bmatrix[0][2 * j] * X[d->elements[i][j] - 1][0] + bmatrix[0][2 * j + 1] * X[d->elements[i][j] - 1][1];
            BX[1] = bmatrix[1][2 * j] * X[d->elements[i][j] - 1][0] + bmatrix[1][2 * j + 1] * X[d->elements[i][j] - 1][1];
            BX[2] = bmatrix[2][2 * j] * X[d->elements[i][j] - 1][0] + bmatrix[2][2 * j + 1] * X[d->elements[i][j] - 1][1];
            stress_x[d->elements[i][j] - 1] = D[0][0] * BX[0] + D[0][1] * BX[1] + D[0][2] * BX[2];
            stress_y[d->elements[i][j] - 1] = D[1][0] * BX[0] + D[1][1] * BX[1] + D[1][2] * BX[2];
            stress_xy[d->elements[i][j] - 1] = D[2][0] * BX[0] + D[2][1] * BX[1] + D[2][2] * BX[2];
        }
    }
    //写入文件
    vector<string> results;

    results.push_back("ZONE N=" + to_string(d->nodes.size()) + ", E=" + to_string(d->elements.size()) + ", F=FEPOINT, ET=QUADRILATERAL");

    for (size_t i = 0; i < d->nodes.size(); ++i)
    {
        results.push_back(to_string(d->nodes[i][0]) + "  " + to_string(d->nodes[i][1]) + "  " + to_string(0) + "  " + to_string(X[i][0]) + "   " + to_string(X[i][1]) + "  " + to_string(stress_x[i]) + "   " + to_string(stress_y[i]) + "  " + to_string(stress_xy[i]));
    }
    for (size_t i = 0; i < d->elements.size(); ++i)
    {
        results.push_back(to_string(d->elements[i][0]) + "  " + to_string(d->elements[i][1]) + "  " + to_string(d->elements[i][2]) + "  " + to_string(d->elements[i][3]));
    }
    //没有输出文件路径时打印到屏幕
    if (argc == 2)
    {
        for (auto _ : results)
        {
            cout << _ << endl;
        }
    }
    else
    {
        writeFile(argv[2], results);
    }
    delete d;
    return 0;
}
