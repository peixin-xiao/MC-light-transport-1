// model generator.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
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

#define Nb 25;  /*Nbins*/
#define N int(512)

struct Tissue {
    int    index;
    char   name[10];
    float miu_a;
    float miu_s;
    float g;
    float index_of_refraction;
};

short int* create3Dmodel(int Nbins) {

    /* USER CHOICES %%%%%%%% <-------- You must set these parameters ------*/
    int SAVEON = 1;        //% 1 = save myname_T.bin, myname_H.mci
    // % 0 = don't save. Just check the program.

    char myname[11] = "skinvessel";//% name for files: myname_T.bin, myname_H.mci
    int time_min = 10;      	//% time duration of the simulation [min] <----- run time -----
    double nm = 532;   	//% desired wavelength of simulation
//    int Nbins       = 20;    	//%  of bins in each dimension of cube
    double binsize = 0.0005; 	//% size of each bin, eg. [cm] or [mm]
//    double T[Nbins][Nbins][Nbins][3];
    short int* T;
    int i, j, k, flag, f = 0;
    int A = 0, B = 0;
    /*sphere center*/

    //    double sphere_center[3] = {12.5,12.5,12.5};
    //    int radius = 5;
    double sphere_center[3] = { Nbins * 0.5,Nbins * 0.5,Nbins * 0.5 };
    int radius = int(Nbins * 0.25);

    T = (short int*)malloc(1 * pow(Nbins, 3) * sizeof(short int));


    printf("success00\n");
    for (i = 0; i < Nbins;i++) {
        for (j = 0;j < Nbins;j++) {
            for (k = 0;k < Nbins;k++) {

                flag = f % 3;
                if (pow(i - sphere_center[0], 2) + pow(j - sphere_center[1], 2) + pow(k - sphere_center[2], 2) >
                    pow(radius, 2)) {
                    *(T + i * Nbins * Nbins + j * Nbins + k) = 0;

                    A++;

                }
                else {
                    *(T + i * Nbins * Nbins + j * Nbins + k) = 1;

                    B++;
                }

                //



            }
            //            printf("success02\n");
        }
        //printf("success03\n");
    }

    //    printf("success04\n");
    //    i = 0; j = Nbins*0.5; k = Nbins*0.5; l = Nbins*0.5;
    printf("%d,%d\n", A, B);
    //    printf("%f,%f,%f\n", *(T + i * Nbins * Nbins * Nbins + j * Nbins * Nbins + k * Nbins + l-1));

    return (T);
}

void write_model(short int* T, std::string directory, int Nbins, float length_voxel) {
    int i, j, k;
    //int Nbins = Nb;
    std::ofstream fout(directory, std::ios::binary);
    fout.write((char*)&(Nbins), sizeof(int));
    fout.write((char*)&(length_voxel), sizeof(float));
    for (i = 0; i < Nbins;i++) {
        for (j = 0;j < Nbins;j++) {
            for (k = 0;k < Nbins;k++) {

                fout.write((char*)&(*(T + i * Nbins * Nbins +
                    j * Nbins + k)), sizeof(short int));

            }
        }
    }
    fout.close();
    std::cout << "model written\n";
}

short int* read_model() {
    int i, j, k;
    int Nbins;
    short int* T;
    std::ifstream fin("model.dat", std::ios::binary);
    fin.read((char*)&Nbins, sizeof(int));
    T = (short int*)malloc(1 * pow(Nbins, 3) * sizeof(short int));
    for (i = 0; i < Nbins;i++) {
        for (j = 0;j < Nbins;j++) {
            for (k = 0;k < Nbins;k++) {

                fin.read((char*)&(*(T + i * Nbins * Nbins +
                    j * Nbins + k)), sizeof(short int));

            }
        }
    }
    fin.close();
    return (T);
    std::cout << "model written\n";
}

std::tuple<int, int(*)[3] > find_vertex(int nbins, short int* t) {
    int i, j, k;
    int vertex[N][3], flag = 0;
    int(*p)[3];
    for (i = 0;i < nbins;i++) {
        for (j = 0;j < nbins;j++) {
            for (k = 0;k < nbins - 1;k++) {
                if (*(t + i * nbins * nbins + j * nbins + k) != *(t + i * nbins * nbins + j * nbins + k + 1))
                {
                    if (std::min_element(t + i * nbins * nbins + j * nbins + k,
                        t + i * nbins * nbins + j * nbins + k + 1) == t + i * nbins * nbins + j * nbins + k + 1)
                    {
                        vertex[flag][0] = i;
                        vertex[flag][1] = j;
                        vertex[flag][2] = k;
                    }
                    else
                    {
                        vertex[flag][0] = i;
                        vertex[flag][1] = j;
                        vertex[flag][2] = k + 1;

                    }
                    flag++;
                }



            }
        }
    }
    auto FLAG = flag;
    p = vertex;
    return std::make_tuple(FLAG, p);
}

void plot_geo(int FLAG, int vertex[N][3]) {
    using namespace matplot;
    std::vector<double> x, y, z;
    int flag;
    x.reserve(FLAG * sizeof(double));
    y.reserve(FLAG * sizeof(double));
    z.reserve(FLAG * sizeof(double));
    for (flag = 0;flag < FLAG;flag++) {
        x.push_back(vertex[flag][0]);
        y.push_back(vertex[flag][1]);
        z.push_back(vertex[flag][2]);
    }

    scatter3(x, y, z);
    show();

}

int main()
{
    int Nbins = Nb;
    float length_voxel = 0.002;
    char directory[] = "C:\\Users\\Administrator\\source\\data\\model.dat";
    short int* T = create3Dmodel(Nbins);
    write_model(T, directory, Nbins, length_voxel);
    int FLAG;
    int(*p)[3];
    std::tie(FLAG, p) = find_vertex(Nbins, T);
    plot_geo(FLAG, p);
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
