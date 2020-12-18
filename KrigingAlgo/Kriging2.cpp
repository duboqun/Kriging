#include"Stdafx.h"
#include"omp.h"
//Standard Output
#include <stdio.h>
#include"Kriging2.h"
#include"GsTL/kriging.h"
#include<GsTL/kriging/OK_constraints.h>


using namespace std;
using namespace System::Runtime::InteropServices;

class Location2d {
public:
	typedef double coordinate_type;

	Location2d() { X = 0;Y = 0; };
	Location2d(int x, int y) : X(x), Y(y) {};
	int& operator[](int i) {
		if (i == 1) return X;
		else     return Y;
	}

	const int& operator[](int i) const {
		if (i == 1) return X;
		else     return Y;
	}

	int X;
	int Y;
};

class geo_value2d {
public:
	typedef double property_type;
	typedef Location2d location_type;

	geo_value2d() :pval_(0) { loc_.X = 0; loc_.Y = 0; };
	geo_value2d(Location2d& loc, double pval) { loc_ = loc; pval_ = pval; };
	geo_value2d(int x, int y, double pval) { loc_.X = x; loc_.Y = y;pval_ = pval; };
	inline Location2d location() const { return loc_; };
	inline double& property_value() { return pval_; }
	inline const double& property_value() const { return pval_; }

private:
	Location2d loc_;
	double pval_;
};

std::ostream& operator << (std::ostream& os, geo_value2d& N) {
	std::cout << N.location().X << " "
		<< N.location().Y << " ";

	return os;
};

typedef std::vector<geo_value2d> Neighborhood;

// gaussian covariance of range 4
inline double gauss_covariance(Location2d u1, Location2d u2) {
	double square_dist = pow(u1[0] - u2[0], 2) +
		pow(u1[1] - u2[1], 2);
	return exp(-3 * square_dist / 16);
}


int maina()
{
	Location2d u1(2, 3);
	Location2d u2(4, -7);

	geo_value2d Z1(u1, 0.21);
	geo_value2d Z2(u2, 0.09);

	Neighborhood neighbors;
	neighbors.push_back(Z1);
	neighbors.push_back(Z2);
	
	Location2d u(0, 0);

	std::vector<double> weights;
	OK_constraints ok;
	double kvariance;
	kriging_weights(weights, kvariance,
		u, neighbors,
		gauss_covariance, ok);
	return 1;
}

/**
 * @brief Calculate the kriging interpolation for a given block of known points given a target resolution
 * @param Z List of known points, the samples, each with the form {x,y,z}
 * @param N size of Z
 * @param RX target horizontal resolution of the interpolated surface
 * @param RY target vertical resolution of the interpolated surface
 * @param VAR semivariogram model to apply
 * @return a Surface with the interpolated points in the target resolution
 */
Surface kriging(vector3 Z[], int N, double RX, double RY, VariogramModel VAR)
{
	//start = clock();
	//initTime();
	Surface S;

	//Get BOUNDS of data minX, maxX, minY, maxY
	double minX = 1000000000, maxX = -1000000000, minY = 1000000000, maxY = -1000000000;
	for (int i = 0; i < N; ++i)
	{
		if (Z[i].x < minX) minX = Z[i].x;
		if (Z[i].x > maxX) maxX = Z[i].x;
		if (Z[i].y < minY) minY = Z[i].y;
		if (Z[i].y > maxY) maxY = Z[i].y;
	}
	S.minX = minX;
	S.minY = minY;
	S.maxX = maxX;
	S.maxY = maxY;

	//Calculate target size in coordinates
	double width = maxX - minX;
	double height = maxY - minY;
	//Calculate target size in points given the spatial resolutions RX, RY, width and height of the target set
	//The size of the target set must include the boundaries of the set given by the known points, that's why an extra row and column are added
	int sizeX = (int)(floor(width / RX)) + 1;
	int sizeY = (int)(floor(height / RY)) + 1;

	cout << "SURFACE DIMENSIONS IN UNITS: " << sizeX << " X " << sizeY << endl;

	S.K = sizeX * sizeY;

	//Create target array and calculate its coordinates
	S.S = new vector3[S.K];

	for (int x = 0; x < sizeX; ++x)
	{
		double coordX = minX + x * RX;
		for (int y = 0; y < sizeY; ++y)
		{
			double coordY = minY + y * RY;
			S.S[x * sizeY + y].x = coordX;
			S.S[x * sizeY + y].y = coordY;
		}
	}

	printTime("PREPARATION");

	//Calculate covariance from each known point to each other
	matrix V = calculateVariogram(Z, Z, N, N, 0, N, 0, N, VAR);

	printTime("Z VARIOGRAM");

	V.invert();
	printTime("V INVERTED");

	//Calculate S
	int n = S.K;
	int q = n / 100;
	int p = q;
	for (int i = 0; i < S.K; ++i)
	{
		matrix D = matrix(N, 1);
		if (q != 0 && i > p)
		{
			p += q;
			int prog = p / q;
			std::cout << "\r" << prog - 1 << "% completed: ";
			std::cout.flush();
		}
		//Calculate Variogram Model for point i
		calculateVariogram(S.S[i], Z, N, D, VAR);
		matrix W = V.I->multiply(D);
		S.S[i].z = 0;
		for (int n = 0; n < N; ++n)
		{
			S.S[i].z += W(n, 0) * Z[n].z;
			//printf("W: %f Z: %f Zo: %f S: %f\n", W[n*N+k], Z[n][2], W[n*N+k] * Z[n][2], S.S[k][2]);
		}
	}

	printTime("KRIGING");
	return S;
}

////////////////////////////////////////////////////

////////////////////////////////////////////////////

/**
 * @brief Main method of the Kriging program, parses a file and calculates the missing points in the surface described by the target resolution
 * @param argc argument count
 * @param argv argument values, in order: filename [nugget sill range variogram resolutionX resolutionY]
 * @return error code, 0 if no errors ocurred
 */
int main2(int argc, char** argv)
{
	VariogramModel v;
	v.C0 = 0.1;
	v.CX = 8.5;
	v.A = 25;
	v.VAR = GAUSSIAN;

	// Test if any arguments were sent, print usage otherwise
	if (argc < 2)
	{
		printf("USAGE: Kriging xyz_file [nugget sill range variogram]\n"); //TODO add to parameters the target resolution values Rx and Ry
		return 0;
	}
	// Check if nugget is sane and store, exit with error otherwise
	if (argc > 2)
	{
		v.C0 = atof(argv[2]);
		if (v.C0 == 0)
		{
			printf("\nNUGGET CANNOT BE 0\n\n");
			return 1;
		}
	}
	// Save Sill
	if (argc > 3)
	{
		v.CX = atof(argv[3]);
	}
	// Save Range
	if (argc > 4)
	{
		v.A = atof(argv[4]);
	}
	// Save Var
	if (argc > 5)
	{
		v.VAR = atoi(argv[5]);
	}
	//TODO capture the target resolution for the final surface Rx and Ry

	printf("FILENAME: %s\n", argv[1]);

	//Load file from disk given the route of the file in ARGV[1]
	int c = 0;
	int e = 0;
	vector3* Z = parseXYZFile(argv[1], c, e);
	if (e == -1)
	{
		printf("\nINPUT FILE IS INCORRECT\n\n");
		return 1;
	}
	else if (e == -2)
	{
		printf("\nINPUT FILE HAS NO LINES\n\n");
		return 1;
	}

	cout << c << endl;

	//TODO send as parameter the target resolution Rx and Ry instead of 1 1
	Surface S = kriging(Z, c, 1, 1, v);

	stringstream ss;
	ss << "Kriging-" << argv[1] << "-" << v.C0 << "-" << v.CX << "-" << v.A << ".xyz";
	writeToFile(S.S, S.K, ss.str().c_str());
	printTime("FILE: - " + ss.str());

	//	printf("\n");
	//
	//	for(int s = 0; s < S.K; ++s)
	//	{
	//		printf("X: %f Y: %f Z: %f\n", S.S[s][0], S.S[s][1], S.S[s][2]);
	//	}

		//delete2DArray(S.S,S.K);
		//delete2DArray(Z,c);

	return 0;
}

KrigingAlgo::Kriging2::Kriging2(cli::array<double>^ xArray, cli::array<double>^ yArray, cli::array<double>^ zArray, RasterContext^ context)
{
	/*
	#define SPHERICAL 0
	#define EXPONENTIAL 1
	#define GAUSSIAN 2
	#define WAVE 3
	#define RATIONAL_Q 4
	#define CIRCULAR 5
	*/
	this->rasterContext = context;
	vecArray = new vector3[xArray->Length];
	vecLength = xArray->Length;
	for (size_t i = 0; i < xArray->Length; i++)
	{
		vecArray[i] = vector3(xArray[i], yArray[i], zArray[i]);
	}
}

KrigingAlgo::Kriging2::~Kriging2()
{
	delete vecArray;
}

void KrigingAlgo::Kriging2::SimpleKrige(Model model, cli::array<float>^ buffer, double nugget, double sill, double range)
{
	VariogramModel v;
	v.C0 = nugget;
	v.CX = sill;
	v.A = range;
	switch (model)
	{
	case KrigingAlgo::Model::Linear:
		v.VAR = CIRCULAR;
		break;
	case KrigingAlgo::Model::LinearWithoutIntercept:
		v.VAR = RATIONAL_Q;
		break;
	case KrigingAlgo::Model::Spherical:
		v.VAR = SPHERICAL;
		break;
	case KrigingAlgo::Model::Exponential:
		v.VAR = EXPONENTIAL;
		break;
	case KrigingAlgo::Model::Gaussian:
		v.VAR = GAUSSIAN;
		break;
	default:
		v.VAR = GAUSSIAN;
		break;
	}

	Surface S;

	S.minX = rasterContext->xOffsetSource;
	S.minY = rasterContext->yOffsetSource;
	S.maxX = rasterContext->xOffsetSource + rasterContext->widthSource * rasterContext->xScaleSource;
	S.maxY = rasterContext->yOffsetSource + rasterContext->heightSource * rasterContext->yScaleSource;

	double RX = rasterContext->xScaleSource;
	double RY = rasterContext->yScaleSource;

	double width = S.maxX - S.minX;
	double height = S.maxY - S.minY;

	int sizeX = rasterContext->widthSource;
	int sizeY = rasterContext->heightSource;


	S.K = sizeX * sizeY;

	S.S = new vector3[S.K];

	for (int x = 0; x < sizeX; ++x)
	{
		double coordX = S.minX + x * RX;
		for (int y = 0; y < sizeY; ++y)
		{
			double coordY = S.minY + y * RY;
			S.S[x * sizeY + y].x = coordX;
			S.S[x * sizeY + y].y = coordY;
		}
	}

	matrix V = calculateVariogram(vecArray, vecArray, vecLength, vecLength, 0, vecLength, 0, vecLength, v);

	V.invert();

	int n = S.K;
	int q = n / 100;
	int p = q;

	GCHandle handle = GCHandle::Alloc(buffer, GCHandleType::Pinned);
	IntPtr bufferPtr = handle.AddrOfPinnedObject();
	float* ptr = static_cast<float*>(bufferPtr.ToPointer());
#pragma omp parallel for
	for (int i = 0; i < S.K; ++i)
	{
		matrix D = matrix(vecLength, 1);
		if (q != 0 && i > p)
		{
			p += q;
			int prog = p / q;
		}
		//Calculate Variogram Model for point i
		calculateVariogram(S.S[i], vecArray, vecLength, D, v);
		matrix W = V.I->multiply(D);
		S.S[i].z = 0;
		for (int n = 0; n < vecLength; ++n)
		{
			S.S[i].z += W(n, 0) * vecArray[n].z;
		}
		ptr[i] = S.S[i].z;
	}
	handle.Free();
}
