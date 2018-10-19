#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NMAX 100
#define MAXTAU 1000
#define PI 3.14

int  fccposn(double x[], double y[], double z[], double alpha, double beta, int Nshell);

int  forces(double x[], double y[], double z[], int N, double alpha, double beta, double L, double V0,
				double Fxi[], double Fyi[], double Fzi[], double Fxe[], double Fye[], double Fze[],
            double *Vint, double *Vext, double rij[][NMAX] );

int  intdyn(double  x[], double  y[], double  z[], double vx[], double vy[], double vz[],
            int N, int nstep, double dt, double gamma, double alpha, double beta,
            double Fxi[], double Fyi[], double Fzi[],  double L, double V0,
            double Fxe[], double Fye[], double Fze[],  double *Vint, double *Vext,
            double rij[][NMAX] );

double ranveladd(double vx[], double vy[], double vz[], int N, double dE, long int *iseed);

double fft1d_(double A[], double work[], int *N, int *isign, int *MNB, double trig[]);

int main(){
	FILE *fcc_file, *annealed_file, *high_energy_file, *cooled_file, *positions, *corr_file, *rij_file, *kinetic_file, *spec;

	double x[NMAX], y[NMAX], z[NMAX]; // arrays holding the positions, velocities and internal and external forces.
	double vx[NMAX], vy[NMAX], vz[NMAX];
	double Fxi[NMAX], Fyi[NMAX], Fzi[NMAX];
	double Fxe[NMAX], Fye[NMAX], Fze[NMAX];
	double alpha, beta; // Lennard Jones parameters
	double L, V0, gamma; // box width, box potential and equations of motion decay parameter.
	double T, nstep, dt; // total time, number of steps and time step
	int N,  Nshell; // number of atoms, number of shells
	int i;
	double *Vint, *Vext; // internal and external potential
	double rij[NMAX][NMAX]; // interatomic distances
	double set, set2; // dummies for initializing pointers
	long int set3, set4, *iseed;
	double E, dE, KE, dKE, E_max; // energy parameters
	int j,h, itau, iav;

	double C_v[MAXTAU], vxold[NMAX][MAXTAU], vyold[NMAX][MAXTAU], vzold[NMAX][MAXTAU]; // velocity correlation arrays
	double A[4*MAXTAU+2], work[4*MAXTAU+2], trig[4*MAXTAU+2]; // fft arrays
	int Nav, maxtau, isign;


	/*open files*/
	fcc_file = fopen("fcc.txt","w"); 				// initial positions of atoms
	annealed_file = fopen("int.txt","w");		// positions of atoms after annealing
	high_energy_file = fopen("mid.txt","w");// positions of atoms after heating
	cooled_file = fopen("fin.txt","w");			// positions of atoms after cooling
	positions = fopen("positions.txt", "w"); // redundant
	corr_file = fopen("corr_file.txt","w");	// correlation file
	rij_file = fopen("rij.txt","w");				// rif file
	kinetic_file = fopen("KE.txt","w");			// kinetic energy during heat cap
	spec = fopen("spectrum.txt", "w");			// fft spectrum

	alpha = 1.0;
	beta = 1.0;
	L = 20.0; // make big enough

	Nshell = 5; /5 shells

	V0 = 0.1;
	vx[0] = 0;
	vy[0] = 0;
	vz[0] = 0;

	gamma = 0;
	nstep = 30000.0;
	dt = 25.0/nstep;

	set = 1.0;
	set2 = 2.0;

	Vint = &set;
	Vext = &set2;

	set3 = 234189;
	set4 = 0;
	iseed = &set3;
	dran1_(&iseed);
	iseed = &set4;

	E = 0.0;
	dKE = 0.0;
	E_max = 3;
	dE = E_max/nstep;

	//=====================================================================================================================
	// Set up the lattice
	N = fccposn(x, y, z, alpha, beta, Nshell); // this adds the atomic positions to x, y, z.
	printf("\nThe number of atoms in the fcc fragment is %d\n\n", N);

	fprintf(positions, "%d\n", N);
	for(i = 0 ; i < N ; i++){
			fprintf(fcc_file, "%lf\t%lf\t%lf\n", x[i], y[i], z[i]);
	}

	//=====================================================================================================================
	//Anneal
	dE = 1;
	dKE = ranveladd(vx, vy, vz, N, dE, iseed); // add energy dE
	intdyn(x, y, z, vx, vy, vz, N, nstep, dt, gamma, alpha, beta, Fxi, Fyi, Fzi, L, V0, Fxe, Fye, Fze, Vint, Vext, rij);

	for(i = 0 ; i < N ; i++){ // print positions after annealing to file
		fprintf(annealed_file, "%lf \t %lf \t %lf \t \n", x[i], y[i], z[i]);
		fprintf(rij_file, "%d \t%lf\t\n", i, rij[i][0]);
	}

	//=====================================================================================================================
	// Heat capacity vs applied energy
	int l = 0;

	// heating up stage.
	for(h=0; h<1000; h++){

		E = E + dE;
		dKE = ranveladd(vx, vy, vz, N, dE, iseed); // add energy then integrate
		intdyn(x, y, z, vx, vy, vz, N, nstep, dt, gamma, alpha, beta, Fxi, Fyi, Fzi, L, V0, Fxe, Fye, Fze, Vint, Vext, rij);

		// calculate kinetic energy
		KE = 0.0;
		for(j=0; j<N; j++){
			KE += (0.5*(vx[j]*vx[j] + vy[j]*vy[j] + vz[j]*vz[j]))/N;
		}
		fprintf(kinetic_file, "%lf\t %lf\n", E , KE);

	}

	for(i = 0 ; i < N ; i++){ // save high energy positions.
		fprintf(high_energy_file, "%lf \t %lf \t %lf \t \n", x[i], y[i], z[i]);
	}

	// cool down stage
	gamma = 0.1;
	for(h=0; h<1000; h++){

		//no energy adding, then integrate
		intdyn(x, y, z, vx, vy, vz, N, nstep, dt, gamma, alpha, beta, Fxi, Fyi, Fzi, L, V0, Fxe, Fye, Fze, Vint, Vext, rij);

		//calculate kinetic energy.
		KE = 0.0;
		for(j=0; j<N; j++){
			KE += (0.5*(vx[j]*vx[j] + vy[j]*vy[j] + vz[j]*vz[j]))/N;
		}
		fprintf(kinetic_file, "%lf\t %lf\n", E , KE);

	}

	for(i = 0 ; i < N ; i++){ // save cooled down positions
		fprintf(cooled_file, "%lf \t %lf \t %lf \t \n", x[i], y[i], z[i]);
	}

	//=====================================================================================================================
	// Vibrational spectrum, the heat cap section should be commented out when this is run.

 	Nav = 5000;
	maxtau = 512;
	dt = 0.02;
	nstep = 5;
	dE = 2;
	gamma = 0;

	dKE = ranveladd(vx, vy, vz, N, dE, iseed); // add a little bit of energy to make harmonic.

	// generate old velocities
	for(i = maxtau-1; i >= 0; i--){

		intdyn(x, y, z, vx, vy, vz, N, nstep, dt, gamma, alpha, beta, Fxi, Fyi, Fzi, L, V0, Fxe, Fye, Fze, Vint, Vext, rij);
		C_v[i] = 0.0;
		for(j = 0; j<N; j++){
			vxold[j][i] = vx[j];
			vyold[j][i] = vy[j];
			vzold[j][i] = vz[j];
		}

	}

	i = 0;
	for(iav = 0; iav < Nav; iav++){

		// integrate
		intdyn(x, y, z, vx, vy, vz, N, nstep, dt, gamma, alpha, beta, Fxi, Fyi, Fzi, L, V0, Fxe, Fye, Fze, Vint, Vext, rij);


		for(itau = maxtau-1; itau>0; itau--){
			for(i=0;i<N;i++){
				// shift old velocities forward
				vxold[i][itau]	= vxold[i][itau-1];
				vyold[i][itau]	= vyold[i][itau-1];
				vzold[i][itau]	= vzold[i][itau-1];
				C_v[itau] += vx[i]*vxold[i][itau] + vy[i]*vyold[i][itau] + vz[i]*vzold[i][itau]; // old vels times new vels

			}

	    for(i=0; i < N; i++){
      	C_v[0] +=   vx[i]*vxold[i][itau] + vy[i]*vyold[i][itau] + vz[i]*vzold[i][itau] ;
	      vxold[i][0] = vx[i];
	      vyold[i][0] = vy[i];
	      vzold[i][0] = vz[i];
			}

	}

	for(i = 0; i < maxtau ; i++){
		C_v[i] = C_v[i]/(N*Nav);  // positional normalization
		fprintf(corr, "%f\n", C_v[i]);
	}

	for(i=0; i<maxtau; i++){ // set up array for fft, need only real parts, fft takes every second place as imaginary
		A[4*maxtau - 2*i] = (A[2*i] = C_v[i]); //real part, C_v[i]
		A[4*maxtau - 2*i+1] = (A[2*i+1] = 0.0); // imag part, 0.
	}

	isign = 1;
	T = maxtau*nstep*dt;
	maxtau = 2*maxtau;
	fft1d_(A, work, &maxtau, &isign, &maxtau, trig); // calculate fft

	for(i=0; i < 2*MAXTAU+1; i++){
		fprintf(spec, "%f, %f\n", 2*i*(2*PI)/T, pow(10000*A[2*i],2)); //square fft and print for vibrational spectrum
	}


	printf("T: %f\n", T);
	printf("Freeing.\n\n");
	fclose(fcc_file);
	fclose(annealed_file);
	fclose(rij_file);
	fclose(kinetic_file);
	fclose(corr);
	fclose(spec);

	return 0;
}

