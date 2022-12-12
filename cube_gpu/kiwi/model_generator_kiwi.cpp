// model_generator_cube.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
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

#define Nb 300;  /*Nbins*/
#define N int(4096) //nums of vertex

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
    double binsize = 0.25e-3; 	//% size of each bin, eg. [cm] or [mm]
//    double T[Nbins][Nbins][Nbins][3];
    short int* T;
    int i, j, k, f = 0;
    int A = 0, B = 0,C=0;
    /*sphere center*/

    //    double sphere_center[3] = {12.5,12.5,12.5};
    //    int radius = 5;
    double sphere_center[3] = { Nbins * 0.5,Nbins * 0.5,Nbins * 0.5 };
    int radius = int(Nbins * 0.25);

    T = (short int*)malloc(1 * pow(Nbins, 3) * sizeof(short int));
    
    //memset(T, 1, sizeof(short)* pow(Nbins, 3));

    printf("success00\n");
    for (i = 0; i < Nbins;i++) {
        for (j = 0;j < Nbins;j++) {
            for (k = 0;k < Nbins;k++) {

                
                if ((pow((i-sphere_center[0]),2)/pow(96,2)+ pow((j - sphere_center[1]), 2) / pow(137, 2)
                    + pow((k - sphere_center[2]), 2) / pow(96, 2)<=1)&& (pow((i - sphere_center[0]), 2) / pow(94, 2)
                        + pow((j - sphere_center[1]), 2) / pow(135, 2)
                        + pow((k - sphere_center[2]), 2) / pow(94, 2) > 1)) {
                    *(T + k * Nbins * Nbins + j * Nbins + i) = 2;

                    A++;

                }
                else if ((pow((i - sphere_center[0]), 2) / pow(94, 2) + pow((j - sphere_center[1]), 2) / pow(135, 2)
                    + pow((k - sphere_center[2]), 2) / pow(94, 2) <= 1)&& (((j>=38)&&(j<=189)&&
                        (pow((i - sphere_center[0]), 2)+ pow((k - sphere_center[2]), 2)<1351))==0)) {
                    *(T + k * Nbins * Nbins + j * Nbins + i) = 3;

                    B++;
                }
                else if((((j >= 38) && (j <= 189) &&
                    (pow((i - sphere_center[0]), 2) + pow((k - sphere_center[2]), 2) < 1351)) == 1))
                {
                    *(T + k * Nbins * Nbins + j * Nbins + i) = 4;
                    C++;
                }
                //



            }
            //            printf("success02\n");
        }
        //printf("success03\n");
    }

    //    printf("success04\n");
    //    i = 0; j = Nbins*0.5; k = Nbins*0.5; l = Nbins*0.5;
    printf("%d,%d,%d\n", A, B,C);
    //    printf("%f,%f,%f\n", *(T + i * Nbins * Nbins * Nbins + j * Nbins * Nbins + k * Nbins + l-1));

    return (T);
}

short int* create3Dmodel_lr(int Nbins,int mag) { //lr: low resolution

    /* USER CHOICES %%%%%%%% <-------- You must set these parameters ------*/
    int SAVEON = 1;        //% 1 = save myname_T.bin, myname_H.mci
    // % 0 = don't save. Just check the program.

    char myname[11] = "skinvessel";//% name for files: myname_T.bin, myname_H.mci
    int time_min = 10;      	//% time duration of the simulation [min] <----- run time -----
    double nm = 532;   	//% desired wavelength of simulation
//    int Nbins       = 20;    	//%  of bins in each dimension of cube
    double binsize = 0.25e-3*mag; 	//% size of each bin, eg. [cm] or [mm]
//    double T[Nbins][Nbins][Nbins][3];
    short int* T;
    int i, j, k, f = 0;
    int A = 0, B = 0, C = 0;
    /*sphere center*/

    //    double sphere_center[3] = {12.5,12.5,12.5};
    //    int radius = 5;
    double sphere_center[3] = { Nbins * 0.5,Nbins * 0.5,Nbins * 0.5 };
    int radius = int(Nbins * 0.25);

    T = (short int*)malloc(1 * pow(Nbins, 3) * sizeof(short int));
    //memset(T, 1, sizeof(short)* pow(Nbins, 3));

    printf("success00\n");
    for (i = 0; i < Nbins;i++) {
        for (j = 0;j < Nbins;j++) {
            for (k = 0;k < Nbins;k++) {


                if ((pow((i - sphere_center[0]), 2) / pow(96/mag, 2) + pow((j - sphere_center[1]), 2) / pow(137/mag, 2)
                    + pow((k - sphere_center[2]), 2) / pow(96/mag, 2) <= 1) && (pow((i - sphere_center[0]), 2) / pow(94/mag, 2)
                        + pow((j - sphere_center[1]), 2) / pow(135/mag, 2)
                        + pow((k - sphere_center[2]), 2) / pow(94/mag, 2) > 1)) {
                    *(T + k * Nbins * Nbins + j * Nbins + i) = 2;

                    A++;

                }
                else if ((pow((i - sphere_center[0]), 2) / pow(94/mag, 2) + pow((j - sphere_center[1]), 2) / pow(135/mag, 2)
                    + pow((k - sphere_center[2]), 2) / pow(94/mag, 2) <= 1) && (((j >= 38/mag) && (j <= 189/mag) &&
                        (pow((i - sphere_center[0]), 2) + pow((k - sphere_center[2]), 2) < pow(36.76/mag,2))) == 0)) {
                    *(T + k * Nbins * Nbins + j * Nbins + i) = 3;

                    B++;
                }
                else if ((((j >= 38/mag) && (j <= 189/mag) &&
                    (pow((i - sphere_center[0]), 2) + pow((k - sphere_center[2]), 2) < pow(36.76 / mag, 2))) == 1))
                {
                    *(T + k * Nbins * Nbins + j * Nbins + i) = 4;
                    C++;
                }
                //



            }
            //            printf("success02\n");
        }
        //printf("success03\n");
    }

    //    printf("success04\n");
    //    i = 0; j = Nbins*0.5; k = Nbins*0.5; l = Nbins*0.5;
    printf("%d,%d,%d\n", A, B, C);
    //    printf("%f,%f,%f\n", *(T + i * Nbins * Nbins * Nbins + j * Nbins * Nbins + k * Nbins + l-1));

    return (T);
}

void write_model(short int* T, std::string directory, int Nbins,
    float length_voxel, int FLAG, int(*p)[3]) {
    int i, j, k;
    int flag;
    //int Nbins = Nb;
    std::ofstream fout(directory, std::ios::binary);
    fout.write((char*)&(Nbins), sizeof(int));
    fout.write((char*)&(length_voxel), sizeof(float));
    for (i = 0; i < Nbins;i++) {
        for (j = 0;j < Nbins;j++) {
            for (k = 0;k < Nbins;k++) {

                fout.write((char*)&(*(T + k * Nbins * Nbins +
                    j * Nbins + i)), sizeof(short int));

            }
        }
    }
    fout.write((char*)&(FLAG), sizeof(int));
    for (flag = 0;flag < FLAG;flag++) {
        fout.write((char*)&(p[flag][0]), sizeof(int));
        fout.write((char*)&(p[flag][1]), sizeof(int));
        fout.write((char*)&(p[flag][2]), sizeof(int));
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
    static int vertex[N][3], flag = 0;
    int(*p)[3];
    for (i = 0;i < nbins;i++) {
        for (j = 0;j < nbins;j++) {
            for (k = 0;k < nbins - 1;k++) {
                if (*(t + k * nbins * nbins + j * nbins + i) != *(t + k * nbins * nbins + j * nbins + i + 1))
                {
                    if (std::min_element(t + k * nbins * nbins + j * nbins + i,
                        t + k * nbins * nbins + j * nbins + i + 1) == t + k * nbins * nbins + j * nbins + i + 1)
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

    p = vertex;
    return std::make_tuple(flag, p);
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
    float length_voxel = 0.25e-3;
    int FLAG;
    int(*p)[3];
    char directory[] = "C:\\Users\\Administrator\\source\\data\\model.dat";
    short int* T = create3Dmodel(Nbins);
    short int* T_lr = create3Dmodel_lr(int(Nbins / 10),10);
    std::tie(FLAG, p) = find_vertex(int(Nbins / 10), T_lr);
    write_model(T, directory, Nbins, length_voxel, FLAG, p);
    plot_geo(FLAG, p);
}

