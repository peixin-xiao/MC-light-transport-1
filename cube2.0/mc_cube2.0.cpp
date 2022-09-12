// mc_cube2.0.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#pragma comment( linker, "/subsystem:windows /entry:mainCRTStartup" )

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <tuple>
#include <random>
#include <iostream>
#include <fstream>
#include <cstdlib> 
#include <ctime>  

#define M 100//amount of photons
#define N int(4096) // amount of vertex of model geometry

#define COUNT 655
#define K int(1e7)

#define	PI          3.1415926
#define	LIGHTSPEED	2.997925E8 /* in vacuo speed of light [m/s] */
#define ALIVE       1   		/* if photon not yet terminated */
#define DEAD        0    		/* if photon is to be terminated */
#define THRESHOLD   0.01		/* used in roulette */
#define CHANCE      0.1  		/* used in roulette */
#define ONE_MINUS_COSZERO 1.0E-12
/* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */
#define SIGN(x)           ((x)>=0 ? 1:-1)

struct Tissue {
    int    index;
    char   name[10];
    float miu_a;
    float miu_s;
    float g;
    float index_of_refraction;
};

std::tuple<float(*)[2]> bubbleSort(float array[][2], int FLAG) {
    float(*a)[2];
    static float ar[N][2];
    memcpy(ar, array, sizeof(array));
    for (int i = 0; i < FLAG - 1; i++) {
        for (int j = 0; j < FLAG - 1 - i; j++) {
            if (ar[j][0] > ar[j + 1][0]) {
                float temp0 = ar[j + 1][0];
                ar[j + 1][0] = ar[j][0];
                ar[j][0] = temp0;
                float temp1 = ar[j + 1][1];
                ar[j + 1][1] = ar[j][1];
                ar[j][1] = temp1;
            }
        }
        //printf("this is %d:", i + 1);

    }
    a = ar;
    return std::make_tuple(a);
}

std::tuple<float, float, float> refra(float face1[], float face2[]
    , float face3[], float refra1, float refra2, float ux, float uy, float uz, float weight) {
    float d1[3], d2[3];
    int i;
    for (i = 0;i < 3;i++) {
        d1[i] = face1[i] - face2[i];
        d2[i] = face2[i] - face3[i];
    }

    float di[3] = { ux,uy,uz };/*light incident direction*/

    float n1[3] = { d1[1] * d2[2] - d1[2] * d2[1]
        ,d1[2] * d2[0] - d1[0] * d2[2]
        ,d1[0] * d2[1] - d1[1] * d2[0] };//normal line of the face : n1 = d1 × d2
    //----------------------------------------------------
    float tempn1 = pow(n1[0], 2) + pow(n1[1], 2)
        + pow(n1[2], 2);
    tempn1 = pow(tempn1, 0.5);
    n1[0] = n1[0] / tempn1;
    n1[1] = n1[1] / tempn1;
    n1[2] = n1[2] / tempn1;//make the norm of vector n1 to 1
    //-----------------------------------------------------
    if (n1[0] * di[0] + n1[1] * di[1] + n1[2] * di[2] < 0)
    {
        for (i = 0;i < 3;i++) n1[i] = -n1[i]; /*inverse normal to make sure
                                                        normal of the face is in the
                                                        same direction as light direction*/
    }

    float n2[3] = { n1[1] * di[2] - n1[2] * di[1]
        ,n1[2] * di[0] - n1[0] * di[2]
        ,n1[0] * di[1] - n1[1] * di[0] };/*normal line between the light direction and
                                                   the normal line of the face : n2 = n1 × light
                                                   direction*/
                                                   //----------------------------------------------------
    tempn1 = pow(n2[0], 2) + pow(n2[1], 2)
        + pow(n2[2], 2);
    tempn1 = pow(tempn1, 0.5);
    n2[0] = n2[0] / tempn1;
    n2[1] = n2[1] / tempn1;
    n2[2] = n2[2] / tempn1;//make the norm of vector n1 to 1
    //-----------------------------------------------------
    float theta1 = acos(n1[0] * di[0] + n1[1] * di[1] + n1[2] * di[2]);
    float theta2 = asin(refra1 * sin(theta1) / refra2);

    float theta;
    if (theta2 - theta1 > 0) theta = theta2 - theta1;
    else theta = 2 * PI + theta2 - theta1;
    /* vector di rotate by angle theta about vector n2*/
    float vrot[3];
    float temp0[3] = { n2[1] * di[2] - n2[2] * di[1]
        ,n2[2] * di[0] - n2[0] * di[2]
        ,n2[0] * di[1] - n2[1] * di[0] };/* temp = n2 × di*/
    float temp1 = n2[0] * di[0] + n2[1] * di[1] + n2[2] * di[2];
    float temp2[3] = { n2[0] * temp1,n2[1] * temp1, n2[2] * temp1 };
    for (i = 0;i < 3;i++) {
        vrot[i] = di[i] * cos(theta) + temp0[i] * sin(theta) + temp2[i] * (1 - cos(theta));
    }

    return std::make_tuple(vrot[0], vrot[1], vrot[2]);
}

std::tuple<int, float, short int*, int, int(*)[3]> read_model(std::string directory) {
    int i, j, k;
    int flag, FLAG;
    int(*p)[3];
    static int p1[N][3];
    int Nbins;
    float length_voxel;
    short int* T;
    std::ifstream fin(directory, std::ios::binary);
    fin.read((char*)&Nbins, sizeof(int));
    fin.read((char*)&length_voxel, sizeof(float));
    T = (short int*)malloc(1 * pow(Nbins, 3) * sizeof(short int));
    for (i = 0; i < Nbins;i++) {
        for (j = 0;j < Nbins;j++) {
            for (k = 0;k < Nbins;k++) {

                fin.read((char*)&(*(T + i * Nbins * Nbins +
                    j * Nbins + k)), sizeof(short int));

            }
        }
    }
    fin.read((char*)&(FLAG), sizeof(int));
    for (flag = 0;flag < FLAG;flag++) {
        fin.read((char*)&(p1[flag][0]), sizeof(int));
        fin.read((char*)&(p1[flag][1]), sizeof(int));
        fin.read((char*)&(p1[flag][2]), sizeof(int));
    }
    fin.close();
    p = p1;
    return std::make_tuple(Nbins, length_voxel, T, FLAG, p);
    std::cout << "model written\n";
}

float randnum(bool type) {
    static float rnd[COUNT];
    static int n;
    float rndvalue = 0;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    if (type == 0) {
        //std::srand(static_cast<unsigned int>(std::time(nullptr)));
        for (int n = 0; n < COUNT; ++n) {
            // Use dis to transform the random unsigned int generated by gen into a 
            // double in [1, 2). Each call to dis(gen) generates a new random double
             rnd[n] = dis(gen);
        }
        
        n = 0;
        rndvalue = rnd[n];
        
    }
    else {
        rndvalue = rnd[n];
        if (n != COUNT - 1) n++;
        else randnum(0);

    }

    return rndvalue;

}

std::tuple<float, float(*)[5], bool, int, float, int, int>distance_to_interface(float x,
    float y, float z, float ux, float uy, float uz, float w,
    short int* T, float delta_t, int Nbins, float s,
    float length_voxel, Tissue a[2], bool reset) {

    float xr = x, yr = y, zr = z;
    float s1 = 0;
    float weight = w;
    static float co[K][5] = { 0 };
    float(*p)[5];
    int flag = 1;
    int i = 0;
    int if_refraction = 0;
    int if_scatter = 0;
    static int time;
    bool if_beyond;
    short int type0 = T[int(x / length_voxel) * Nbins * Nbins
        + int(y / length_voxel) * Nbins + int(z / length_voxel)];
    if_beyond = 0;
    float refraction = a[type0].index_of_refraction;
    float mua = a[type0].miu_a;
    float step = LIGHTSPEED * delta_t / refraction;
    float length_of_side = Nbins * length_voxel;

    if (reset) {
        time = 0;
        memset(co, 0, sizeof(co));
    }

    do {

        if (((xr / length_voxel) - 1) < 0 || ((yr / length_voxel) - 1) < 0 
            || ((zr / length_voxel) - 1) < 0 || 
            ((xr / length_voxel) - 1) > Nbins-1 ||
            ((yr / length_voxel) - 1) > Nbins-1 ||
            ((zr / length_voxel) - 1) > Nbins-1)
        {
            flag = 0;
            if_beyond = 1;
            std::cout << "out of cube\n";
        }
        else
        {

            if (T[int((xr + step * ux) / length_voxel) * Nbins * Nbins +
                int((yr + step * uy) / length_voxel) * Nbins +
                int((zr + step * uz) / length_voxel)] == type0)
            {
                s1 += step;
                xr = xr + step * ux;
                yr = yr + step * uy;
                zr = zr + step * uz;
                weight = weight * (exp(-mua * step));
                co[i][0] = xr;
                co[i][1] = yr;
                co[i][2] = zr;
                co[i][3] = weight;
                co[i][4] = delta_t * time;
                time++;
                i++;
            }
            else {
                //s1 += length_voxel;
                int la = 1;
                xr = xr + la*step * ux;
                yr = yr + la*step * uy;
                zr = zr + la*step * uz;
                co[i][0] = xr;
                co[i][1] = yr;
                co[i][2] = zr;
                co[i][3] = weight;
                co[i][4] = delta_t * time;
                time++;
                i++;
                flag = 0;
                if_refraction = 1;
                std::cout << "reach the interface\n";

            }
        }
        if (s1 >= s) {
            flag = 0;
            if_scatter = 1;
            std::cout << "out of max length\n";
        }
        //std::cout << "flag = " << flag << "\n";
    } while (flag);
    p = co;
    std::cout << "i = " << i << "\n";
    return std::make_tuple(s1, p, if_beyond, i,
        weight, if_refraction, if_scatter);
}

std::tuple<std::vector< std::vector<float> > > record_path(float(*c0)[5], int k, int reset) {
    //static int k0 = 0;

    static std::vector< std::vector<float> > path(5);
    //static float c00[4];
    if (reset == 0) {
        int i, j;
        for (i = 0;i < path.size();i++) {
            //path[i].resize((k0 + k) * sizeof(float));
            for (j = 0;j < k;j++) {
                path[i].push_back(c0[j][i]);
            }
            //c00[i] = c0[j][i];
        }

    }
    else if (reset == 1) {
        path.clear();
        path.resize(5 * sizeof(float));
        //static std::vector< std::vector<float> > path(4);
        std::cout << path.size() << "\n";

        /*
        int i, j;
        for (i = 0;i < 4;i++) {
            for (j = 0;j < k-1;j++) {
                path[i].push_back(float(0.000));
            }


        }
       */

    }

    return make_tuple(path);

    //k0 = k0 + k;

}



std::tuple< std::vector< std::vector<float> >(*) >  montecarlo
(short int* T, int Nbins, float length_voxel, Tissue a[2], int FLAG, int(*vertex)[3]) {

    /* propagation parameters */
    double	x, y, z;        /* photon position */
    double	ux, uy, uz;     /* photon trajectory as cosines */
    double  uxx, uyy, uzz;	/* temporary values used during spin */
    double	s;              /* step sizes. s = -log(rnd)/mus [cm] */
    double  sleft = 0;          /* dimensionless */
    double	costheta;       /* cos(theta) */
    double  sintheta;       /* sin(theta) */
    double	cospsi;         /* cos(psi) */
    double  sinpsi;         /* sin(psi) */
    double	psi;            /* azimuthal angle */
    long	i_photon;       /* current photon */
    double	w;              /* photon weight */
    double	absorb;         /* weighted deposited in a step due to absorption */
    short   photon_status;  /* flag = alive=1 or dead=0 */

    /* other variables */
    double	mua;            /* absorption coefficient [cm^-1] */
    double	mus;            /* scattering coefficient [cm^-1] */
    double	g;              /* anisotropy [-] */
    double  refraction;     /* index of refraction*/
    double	nphotons;       /* number of photons in simulation */

    /* dummy variables */
    double  rnd;            /* assigned random value 0-1 */

    int i_voxel;
    bool flag1 = 0;
    double 	temp;           /* dummy variable */

    /* mcxyz bin variables */
    //float   length_voxel = 0.002;/*bins size, unit : [m] */
    //float   step = 0.5;
    //float   stepsize = length_voxel * step;
    double   delta_t = 1e-14;

    double s1;
    double s3; //distance from initial point to the current position
    //int s_block;
    int if_refraction = 0;
    int if_scatter = 0;

    float(*c0)[5] = { 0 };
    static std::vector< std::vector<float>>  path1[M];
    std::vector< std::vector<float> >* path;
    std::vector< std::vector<float> >  path0;
    int if_beyond, k;
    int i, j;
    int rndflag = 0;

    float dis[N][2], (*dis2)[2]; // all surface point before and after bubble sort
    float face1[3], face2[3], face3[3]; // coordinates of the 3 nearest point
    int pos; //position of the minimum 3 point near the light penetration point
    float refra0, refra1; //refraction index before and after refraction
    //InitRandomGen;
    randnum(0);
    //std::srand(static_cast<unsigned int>(0.5*32767));
    nphotons = M - 1; // will be updated to achieve desired run time, time_min.
    i_photon = -1;

    do {
        i_photon += 1;				/* increment photon count */
        w = 1.0;                    /* set photon weight to one */
        photon_status = ALIVE;      /* launch an alive photon */

        x = 0.5 * Nbins * length_voxel;
        y = 1 * length_voxel;
        z = 0.5 * Nbins * length_voxel;
        ux = 0;
        uy = 1;
        uz = 0;
        sleft = 0;
        std::cout << "\n\n\nphoton : " << i_photon << " launched\n";
        int loop = 0;
        distance_to_interface(x, y, z, ux, uy, uz, w, T,
            delta_t, Nbins, 0, length_voxel, a, 1);
        //std::cout << "distance to interface initialized\n";
        //record_path(c0, k, 1);

        do {
            loop++;
            //if (loop >= 100) break;
            std::cout << "loop = " << loop << " \n";

            if (sleft == 0 || flag1 == 1) {
                //while ((rnd = RandomNum) <= 0.0);
                rnd = randnum(1);
                //rnd = RND;
                sleft = -log(rnd);
                std::cout << "sleft = " << sleft << "\n";
                flag1 = 0;
            }
            else {
                i_voxel = int(x / length_voxel) * Nbins * Nbins
                    + int(y / length_voxel) * Nbins + int(z / length_voxel);
                mua = a[T[i_voxel]].miu_a;
                mus = a[T[i_voxel]].miu_s;
                g = a[T[i_voxel]].g;
                refraction = a[T[i_voxel]].index_of_refraction;

                s = sleft / (mus);  //[m]
                //s_block = s / length_voxel;  // [unit]

                std::tie(s1, c0, if_beyond, k, w,
                    if_refraction, if_scatter) =
                    distance_to_interface(x, y, z, ux, uy, uz, w, T,
                        delta_t, Nbins, s, length_voxel, a, 0);
                //std::cout << "\n\nw = " << w << "\n\n";

                /*i_voxel = int(round((x / length_voxel) * Nbins * Nbins +
                    (y / length_voxel) * Nbins + (z / length_voxel)));*/
                refra0 = a[T[i_voxel]].index_of_refraction;
                tie(path0) = record_path(c0, k, 0);


                path1[i_photon] = path0;


                if (if_beyond == 1) { photon_status = DEAD; }
                if (s1 > s) { s3 = s; flag1 = 1; }
                else s3 = s1;

                //std::cout << "s3 = " << s3 << " " << "s1 = " << s1 << "\n";

                k--;
                //if (k > 0) k--;
                /*int la = 1;
                x = c0[k][0] + la * length_voxel * ux;
                y = c0[k][1] + la * length_voxel * uy;
                z = c0[k][2] + la * length_voxel * uz;*/
                x = c0[k][0];
                y = c0[k][1];
                z = c0[k][2];

                /*i_voxel = int(floor(((x / length_voxel)-1) * Nbins * Nbins +
                    ((y / length_voxel)-1) * Nbins + ((z / length_voxel)-1)));*/
                i_voxel = int(x / length_voxel) * Nbins * Nbins
                    + int(y / length_voxel) * Nbins + int(z / length_voxel);

                //std::cout << "xyz=" << x << y << z << "\n";
                if(i_voxel>pow(Nbins,3)-1) photon_status = DEAD;
                if (x <= 0 || y <= 0 || z <= 0 || x >=
                    Nbins * length_voxel || y >= Nbins * length_voxel ||
                    z >= Nbins * length_voxel) {
                    photon_status = DEAD;
                }
                if (i_voxel > sizeof(T) / sizeof(T[0])) i_voxel = sizeof(T) / sizeof(T[0]);
                refra1 = a[T[i_voxel]].index_of_refraction;
                //absorb = w * (1 - exp(-mua * s3));   /* photon weight absorbed at this step */
                //w -= absorb;
                sleft = sleft - s3 * mus;
                /* decrement weight by amount absorbed */
                std::cout << "sleft = " << sleft << "\n";
            }

            if (if_refraction == 1) {
                for (i = 0;i < FLAG;i++) {
                    dis[i][0] = pow((pow(x - length_voxel * vertex[i][0], 2)
                        + pow(y - length_voxel * vertex[i][1], 2)
                        + pow(z - length_voxel * vertex[i][2], 2)), 0.5);
                    dis[i][1] = i;
                }
                std::tie(dis2) = bubbleSort(dis, FLAG);
                pos = dis2[N - FLAG][1];
                std::cout << "pos = " << pos << std::endl;
                for (i = 0;i < 3;i++) {
                    face1[i] = length_voxel * vertex[pos][i];
                }
                pos = dis2[N - FLAG + 1][1];
                for (i = 0;i < 3;i++) {
                    face2[i] = length_voxel * vertex[pos][i];
                }
                pos = dis2[N - FLAG + 2][1];
                for (i = 0;i < 3;i++) {
                    face3[i] = length_voxel * vertex[pos][i];
                }
                std::make_tuple(ux, uy, uz) =
                    refra(face1, face2, face3, refra0, refra1, ux, uy, uz, w);

            }

            if (if_scatter == 1) {
                /* sample for costheta */
                //rnd = RandomNum;
                rnd = randnum(1);
                rnd = rnd;
                //rnd = RND;
                if (g == 0.0)
                    costheta = 2.0 * rnd - 1.0;
                else {
                    temp = (1.0 - g * g) / (1.0 - g + 2 * g * rnd);
                    costheta = (1.0 + g * g - temp * temp) / (2.0 * g);
                }
                sintheta = sqrt(1.0 - costheta * costheta); /* sqrt() is faster than sin(). */

                /* sample psi. */
                rnd = randnum(1);
                //rnd = RND;
                psi = 2.0 * PI * rnd;
                //psi = 2.0 * PI * RandomNum;
                cospsi = cos(psi);
                if (psi < PI)
                    sinpsi = sqrt(1.0 - cospsi * cospsi);     /* sqrt() is faster than sin(). */
                else
                    sinpsi = -sqrt(1.0 - cospsi * cospsi);

                /* new trajectory. */
                if (1 - fabs(uz) <= ONE_MINUS_COSZERO) {      /* close to perpendicular. */
                    uxx = sintheta * cospsi;
                    uyy = sintheta * sinpsi;
                    uzz = costheta * SIGN(uz);   /* sign() is faster than division. */
                }
                else {					/* usually use this option */
                    temp = sqrt(1.0 - uz * uz);
                    uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
                    uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
                    uzz = -sintheta * cospsi * temp + uz * costheta;
                }

                /* update trajectory */
                ux = uxx;
                uy = uyy;
                uz = uzz;

                std::cout << "new direction = " << ux << "," << uy << "," << uz << "\n";

            }
            /**** check roulette
                 if photon weight below threshold, then terminate photon using roulette technique.
                 photon has chance probability of having its weight increased by factor of 1/chance,
                 and 1-chance probability of terminating.
                 *****/
            if (w < THRESHOLD) {
                rnd = randnum(1);
                //rnd = RND;
                if (rnd <= CHANCE)
                    w /= CHANCE;
                else
                {
                    photon_status = DEAD;
                    //std::cout << "i_photon = " << i_photon << "\n";
                    //std::cout << "coordinate = " << x << "," << y << "," << z << "\n";

                }
            }






        } while (photon_status == ALIVE);  /* end step_check_hop_spin */
        /* if alive, continue propagating */
        /* if photon dead, then launch new photon. */

        for (j = 0;j < k;j++) {
            for (i = 0;i < 5;i++) {
                std::cout << path0[i][j] << " ";
            }
            std::cout << "\n";
        }

        record_path(c0, k, 1);

    } while (i_photon < nphotons);  /* end run */
    path = path1;

    return std::make_tuple(path);

}

void write_path(std::vector< std::vector<float> >(*path), std::string directory) {
    std::ofstream fout(directory, std::ios::binary);
    int i, j, k;
    int I, J, L;
    I = M;
    J = 5;
    L = 50;
    fout.write((char*)&(I), sizeof(int));
    fout.write((char*)&(J), sizeof(int));
    fout.write((char*)&(L), sizeof(int));
    for (i = 0;i < M;i++) {
        for (j = 0;j < J;j++) {
            for (k = 0;k < path[i][j].size();k++) {
                fout.write((char*)&path[i][j][k], sizeof(float));
                //std::cout << path[i][j][k] << "\n";
            }
            fout.write("end", 4);
        }
    }

    fout.close();
}


int main()
{
    Tissue a[2];
    //a[0] = { 0,"air",0.0001,1,1.0 ,1};
    a[0] = { 0,"air",0.01,0.01,0 ,1 };
    a[1] = { 1,"medium 1",0.1,0.24e4,0.7 ,1.6 };//mua ~= 0.001/cm
    int Nbins;
    float length_voxel;
    int(*vertex)[3];
    int FLAG;
    short int* T1;
    std::string directory0 = "C:\\Users\\Administrator\\source\\data\\model.dat";
    std::string directory1 = "C:\\Users\\Administrator\\source\\data\\path.dat";
    std::vector< std::vector<float> >(*path);
    std::tie(Nbins, length_voxel, T1,
        FLAG, vertex) = read_model(directory0);

    std::tie(path) = montecarlo(T1, Nbins, length_voxel, a, FLAG, vertex);
    write_path(path, directory1);
    return 0;
}
