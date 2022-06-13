// ceshi3.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//


// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单
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

#define ls          1.0e-7;
#define	PI          3.1415926
#define	LIGHTSPEED	2.997925E10 /* in vacuo speed of light [cm/s] */
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
#define M 512
/* DECLARE FUNCTION */
struct Tissue {
    int    index;
    char   name[10];
    double miu_a;
    double miu_s;
    double g;
};

double RandomGen(char Type, long Seed, long* Status);
/* Random number generator */


short int* create3Dmodel(int Nbins);


std::tuple<float (*)[3]> montecarlo(short int* T, int Nbins);

bool SameVoxel(double x1, double y1, double z1, double x2, double y2, double z2, double dx, double dy, double dz);

double min2(double a, double b);

void plot_geo(int FLAG, int vertex[N][3]);

void plot_geo2(int FLAG, float vertex[N][3]);

std::tuple<int, int (*)[3]> find_vertex(int nbins, short int* t);

std::tuple<float, float, int(*)[3], int>find_distance_to_boundry(int x, int y, int z, float u1,
    float u2, float u3, short int* T, float step, int Nbins, float s);





int main() {
    
    using namespace matplot; 
    using namespace std;
    Tissue a[2];
    a[0] = { 0,"air",0.0001,1,1.0 };
    a[1] = { 1,"kiwifruit",10,1e4,0.7 };
    int Nbins = 25;
    int flag = M;
    int i = 0;
    float(*q)[3] = { 0 };
    printf("success04\n");
    float rnd;
    short int* T = create3Dmodel(Nbins);
    std::cout << "3D model created\n";
    tie(q) = montecarlo(T, Nbins);
    plot_geo2(flag, q);
    
    



    /*
    int FLAG;
    int(*p)[3];
    std::tie(FLAG, p) = find_vertex(nbins, t);
    plot_geo(FLAG, p);
    */
    /*
    float s1, s2;
    int (*p)[3];
    int x = 0,y=0,z=0;
    float ux = 100, uy = 80, uz = 50;
    float u;
    float step = 0.5;
    int i,j,n;
    float s = 18;

    u = pow(pow(ux, 2) + pow(uy, 2) + pow(uz, 2), 0.5);
    ux = ux / u;uy = uy / u;uz = uz / u;
    std::tie(s1, s2, p) = find_passed_block(x, y, z, ux, uy, uz, T, step, Nbins,s);
    
    std::cout << s1 << " " << s2 << "\n";
    for (i = 0;i < 64;i++) {
        for (j = 0;j < 3;j++) {
            n = p[i][j];
            std::cout << n << " ";
        }
        std::cout << "\n";
    }
   */
    
    /*
    int b;
    int(*p)[3];
    int k;
    tie(b, p) = test1();
    k = p[0][0];
    std::cout << "k=p[0][0]=" << k << "\n";
    std::cout <<"p[0][:]="<< p[0][0] << " " << p[0][1] << "  " << p[0][2] << "\n";
    k = p[0][0];
    std::cout << "k=p[0][0]=" << k << "\n";
    */
    return 0;
    
    
    

    



}





std::tuple<float,float, int(*)[3], int>find_distance_to_boundry(int x, int y, int z, float u1,
    float u2, float u3, short int* T, float step, int Nbins, float s_block) {
    
    
    float xr=x,yr=y,zr=z;
    float s1=0, s2=0;
    static int co[64][3]={0};
    int(*p)[3];
    int i=0,flag=1;
    int f;
    short int type0 = T[x * Nbins * Nbins + y * Nbins + z];
    f = 0;
    
    
    do {

        /*
        std::cout << "Program is paused !\n" <<
            "Press Enter to continue\n";

        // pause the program until user input
        f = getc(stdin);

        std::cout << "\nContinuing .";
        */

        /*if (s1>=s) {
            flag = 0;
            std::cout << "out of up limitation\n";
        }*/
        if (((xr - x) >= 1) || ((yr - y) >= 1) || ((zr - z) >= 1)) {
            x = int(xr);y = int(yr);z = int(zr);
            co[i][0] = x; co[i][1] = y;co[i][2] = z;
            //std::cout << co[i][0] << " " << co[i][1] << " " << co[i][2] << "\n";
            i++;
            //std::cout << "fill in " << i << "\n";
        }
        else {
            xr = xr + step * u1;
            yr = yr + step * u2;
            zr = zr + step * u3;
            //std::cout << "step forward\n";
            //std::cout << "xr=" << xr << " yr=" << yr << " zr=" << zr << "\n";
            if (xr <= 0 || yr <= 0 || zr <= 0 || xr >= Nbins - 1 ||
                yr >= Nbins - 1 || zr >= Nbins - 1)
            {
                flag = 0;
                f = 1;
                std::cout << "out of range\n";
            }
            else
            {
                if (T[int(x * Nbins * Nbins + y * Nbins + z)] == type0)
                {
                    s1 += step;
                }
                else {
                    s2 += step;
                    flag = 0;
                    std::cout << "reach the boundry\n";
                }
            }
        }
        if (s1 >= s_block) {
            flag = 0;
            std::cout << "out of up limitation\n";
        }
    } while (flag);
    p = co;
    //std::cout << p[0][0] << " " << p[0][1] << " " << p[0][2] << "\n";
    return std::make_tuple(s1, s2, p, f);

}

std::tuple<int, int (*)[3] > find_vertex(int nbins, short int* t) {
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


void plot_geo2(int FLAG, float vertex[M][3]) {
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

std::tuple<float (*)[3]> montecarlo(short int* t, int nbins) {
    using namespace std;
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
    bool sv;             /* are they in the same voxel? */

    /* other variables */
    double	mua;            /* absorption coefficient [cm^-1] */
    double	mus;            /* scattering coefficient [cm^-1] */
    double	g;              /* anisotropy [-] */
    double	nphotons;       /* number of photons in simulation */
    short int type;

    /* launch parameters */
    int		mcflag, launchflag, boundaryflag;
    float	xfocus, yfocus, zfocus;
    float	ux0, uy0, uz0;
    float	radius;
    float	waist;

    /* dummy variables */
    double  rnd;            /* assigned random value 0-1 */
    double	r, phi;			/* dummy values */
    int j;
    int i_voxel;
    bool flag1=0;
    long	nn;         /* dummy indices */
    double	tempx, tempy, tempz; /* temporary variables, used during photon step. */
    int 	ix, iy, iz;     /* added. used to track photons */
    double 	temp;           /* dummy variable */
    int     bflag;          /* boundary flag:  0 = photon inside volume. 1 = outside volume */
    int		cnt;

    /* mcxyz bin variables */
    float	dx, dy, dz;     /* bin size [cm] */
    int		nx, ny, nz, nt; /* # of bins */
    float	xs, ys, zs;		/* launch position */
    float   length_voxel = 0.002;/*bins size, unit : [m] */
    float   step = 0.5;
    float   stepsize = length_voxel * step;
    /* time */
    float	time_min;               // requested time duration of computation.
    time_t	now;
    double	start_time, finish_time, temp_time; /* for clock() */

    struct Tissue a[2];
    a[0] = { 0,"air",0.00001,1,0.9 };
    a[1] = { 1,"kiwifruit",0.1,10,0.7 };

    float s1, s2;
    float s3; //distance from initial point to the current position
    int s_block;
    int(*p)[3];
    static float coordinate[M][3] = {0};
    float(*q)[3];
    int f;
    
    
    InitRandomGen;
    nphotons = M; // will be updated to achieve desired run time, time_min.
    i_photon = 0;
    cnt = 0;
    do {
        i_photon += 1;				/* increment photon count */
        w = 1.0;                    /* set photon weight to one */
        photon_status = ALIVE;      /* launch an alive photon */
        cnt = 0;
        x = 1;
        y = 1 ;
        z = 0.5*nbins ;   
        ux = pow(2,0.5);
        uy = pow(2,0.5);
        uz = 0;
        sleft = 0;
        std::cout << "\n\n\nphoton : " << i_photon << " launched\n";
        int loop = 0;
        
        

        do {
            loop++;
            if (loop >= 100) break;
            std::cout << "loop = " << loop << " \n";

            /*if (x<=0 || y<=0 || z<=0 || x>=nbins-1 || y>=nbins-1 || z>=nbins-1) {
                photon_status = DEAD;
            }*/
            if (sleft < 0.001||flag1==1) {
                while ((rnd = RandomNum) <= 0.0);
                sleft = -log(rnd);
                std::cout << "sleft = " << sleft << "\n";
                flag1 = 0;
            }
            else {
                i_voxel = int(x) * nbins * nbins + int(y) * nbins + int(z);
                mua = a[t[i_voxel]].miu_a;
                mus = a[t[i_voxel]].miu_s;
                g = a[t[i_voxel]].g;


                s = sleft / (mus);  //[m]
                s_block = s / length_voxel;  // [unit]

                std::tie(s1, s2, p, f) = find_distance_to_boundry(x, y, z, ux, uy, uz, t, step, nbins, s_block);
                if (f == 1) { photon_status = DEAD; }
                if (s1 > s_block) { s3 = s_block; flag1 = 1; }
                else s3 = s1;
                //std::cout << "s3 = " << s3 <<" "<<"s1 = "<<s1<<" "<<"s_block = "<<s_block<< "\n";
                x = x + s3 * ux;				/* update positions. [cm] */
                y = y + s3 * uy;
                z = z + s3 * uz;
                if (x<=0 || y<=0 || z<=0 || x>=nbins-1 || y>=nbins-1 || z>=nbins-1) {
                photon_status = DEAD;
            }
                absorb = w * (1 - exp(-mua * s3));   /* photon weight absorbed at this step */
                w -= absorb;
                sleft = sleft - s3 * mus * length_voxel;
                /* decrement weight by amount absorbed */
                std::cout << "sleft = " << sleft << "\n";
            }
            if (flag1==1) {
                /* sample for costheta */
                rnd = RandomNum;
                if (g == 0.0)
                    costheta = 2.0 * rnd - 1.0;
                else {
                    double temp = (1.0 - g * g) / (1.0 - g + 2 * g * rnd);
                    costheta = (1.0 + g * g - temp * temp) / (2.0 * g);
                }
                sintheta = sqrt(1.0 - costheta * costheta); /* sqrt() is faster than sin(). */

                /* sample psi. */
                psi = 2.0 * PI * RandomNum;
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
                if (RandomNum <= CHANCE)
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
        
        coordinate[i_photon - 1][0] = x;
        coordinate[i_photon - 1][1] = y;
        coordinate[i_photon - 1][2] = z;
        std::cout << "coordinate = " << x << "," << y << "," << z << "\n";
        std::wcout << "weight = " << w << "\n";
    } while (i_photon < nphotons);  /* end run */
    q = coordinate;
    return std::make_tuple(q);
    
}

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
    int i, j, k, l, flag, f = 0;
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

double RandomGen(char Type, long Seed, long* Status) {
    static long i1, i2, ma[56];   /* ma[0] is not used. */
    long        mj, mk;
    short       i, ii;

    if (Type == 0) {              /* set seed. */
        mj = MSEED - (Seed < 0 ? -Seed : Seed);
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;
        for (i = 1; i <= 54; i++) {
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ)
                mk += MBIG;
            mj = ma[ii];
        }
        for (ii = 1; ii <= 4; ii++)
            for (i = 1; i <= 55; i++) {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ)
                    ma[i] += MBIG;
            }
        i1 = 0;
        i2 = 31;
    }
    else if (Type == 1) {       /* get a number. */
        if (++i1 == 56)
            i1 = 1;
        if (++i2 == 56)
            i2 = 1;
        mj = ma[i1] - ma[i2];
        if (mj < MZ)
            mj += MBIG;
        ma[i1] = mj;
        return (mj * FAC);
    }
    else if (Type == 2) {       /* get status. */
        for (i = 0; i < 55; i++)
            Status[i] = ma[i + 1];
        Status[55] = i1;
        Status[56] = i2;
    }
    else if (Type == 3) {       /* restore status. */
        for (i = 0; i < 55; i++)
            ma[i + 1] = Status[i];
        i1 = Status[55];
        i2 = Status[56];
    }
    else
        puts("Wrong parameter to RandomGen().");
    return (0);
}

bool SameVoxel(double x1, double y1, double z1, double x2, double y2, double z2, double dx, double dy, double dz)
{
    double xmin = min2((floor)(x1 / dx), (floor)(x2 / dx)) * dx;
    double ymin = min2((floor)(y1 / dy), (floor)(y2 / dy)) * dy;
    double zmin = min2((floor)(z1 / dz), (floor)(z2 / dz)) * dz;
    double xmax = xmin + dx;
    double ymax = ymin + dy;
    double zmax = zmin + dz;
    bool sv = 0;

    sv = (x1 <= xmax && x2 <= xmax && y1 <= ymax && y2 <= ymax && z1 < zmax&& z2 <= zmax);
    return (sv);
}

double max2(double a, double b) {
    double m;
    if (a > b)
        m = a;
    else
        m = b;
    return m;
}

double min2(double a, double b) {
    double m;
    if (a >= b)
        m = b;
    else
        m = a;
    return m;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC