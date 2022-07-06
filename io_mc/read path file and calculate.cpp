// read path file and calculate.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <matplot/matplot.h>
#include <tuple>
#include <random>
#include <iostream>
#include <cstdlib> 
#include <ctime> 

#define ls          1.0e-7;
#define	PI          3.1415926
#define	LIGHTSPEED	2.997925E8 /* in vacuo speed of light [m/s] */
#define ALIVE       1   		/* if photon not yet terminated */
#define DEAD        0    		/* if photon is to be terminated */
#define THRESHOLD   0.01		/* used in roulette */
#define CHANCE      0.1  		/* used in roulette */
#define COS90D      1.0E-6
/* If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad. */
#define ONE_MINUS_COSZERO 1.0E-12
/* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */
#define SIGN(x)           ((x)>=0 ? 1:-1)
#define InitRandomGen    (double) RandomGen(0, 1, NULL)
/* Initializes the seed for the random number generator. */
#define RandomNum        (double) RandomGen(1, 0, NULL)
/* Calls for a random number from the randum number generator. */
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9
#define N int(512)
#define M 64 /*number of photons*/
#define K 128 /*the number of time resolved 4-D coordinate returned by function
                            distance_to_interface()*/
#define Nb 25;  /*Nbins*/
#define COUNT 5000  /*quantity of random number*/
#define RND float(std::rand() )/ RAND_MAX

std::tuple<std::vector< std::vector<float> >(*)> read_path(std::string directory) {
    int i, j, k;
    int I, J, L;
    char f[4];
    char end = 'e';
    char f1[] = "end";
    float p;
    std::ifstream fin(directory, std::ios::binary);
    static std::vector<std::vector<float>>  path[M];
    std::vector<std::vector<float>>* path1;
    fin.read((char*)&I, sizeof(int));
    fin.read((char*)&J, sizeof(int));
    fin.read((char*)&L, sizeof(int));
    std::cout << "I = " << I << "  J = " << J << "  L = " << L << "\n";
    for (i = 0;i < I;i++) {
        path[i].resize(J);
        for (j = 0;j < J;j++) {
            for (k = 0;;k++) {
                fin.read((char*)&f, 4);
                if (strcmp(f, "end") != 0) {
                    fin.seekg(-4L, std::ios::cur);
                    fin.read((char*)&p, sizeof(float));
                    path[i][j].push_back(p);
                }
                else {
                    //fin.seekg(3L, std::ios::cur);
                    break;
                }
            }
        }

    }
    fin.close();
    path1 = path;
    return std::make_tuple(path1);
}

void plot_tragectory2(std::vector< std::vector<float> >* path) {
    using namespace matplot;
    int i, j;
    std::vector<std::vector<double>> x(M), y(M), z(M);
    std::cout << "01\n";
    //int k = path[0][0].size() ;
    std::cout << "02\n";
    for (i = 0;i < M;i++) {
        for (j = 0;j < path[i][0].size();j++) {
            x[i].push_back(path[i][0][j]);
            //std::cout << "\n\n\nx["<<i<<"][" << j << "] = " << path[i][0][j] << "\n";
            y[i].push_back(path[i][1][j]);
            z[i].push_back(path[i][2][j]);
        }

    }

    plot3(x, y, z);


    show();
}

int main()
{
    //read the path file and store it in a vector[M] [5][?]
    std::string directory = "C:\\Users\\Administrator\\source\\data\\path.dat";
   std::vector< std::vector<float> >(*path);
   std::tie(path) = read_path(directory);
   plot_tragectory2(path);
   //--------------------------------------------------------------------
    
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
