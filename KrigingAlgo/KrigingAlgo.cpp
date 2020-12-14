// 这是主 DLL 文件。

#include "stdafx.h"

#include "KrigingAlgo.h"

// Kriging.cpp : Defines the entry point for the console application.
//

#include <algorithm>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <assert.h>
#include <direct.h>
using namespace System::Runtime::InteropServices;

// Various helper routines for debugging

// Output a matrix

std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& matrix)
{
	for (size_t i = 0; i < matrix.size(); i++)
	{
		for (size_t j = 0; j < matrix.size(); j++)
		{
			if (j != 0)
			{
				std::cout << " ";
			}

			std::cout << std::setw(10) << std::setprecision(7) << matrix[i][j];
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;

	return os;
}

// Matrix multiplication

std::vector<std::vector<double>> MatrixTimesMatrix(const std::vector<std::vector<double>>& A,
	const std::vector<std::vector<double>>& B)
{
	std::vector<std::vector<double>> C(A.size(), std::vector<double>(B[0].size(), 0.0));

	for (size_t i = 0; i < A.size(); i++)
	{
		for (size_t j = 0; j < B[0].size(); j++)
		{
			for (size_t k = 0; k < A[0].size(); k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return C;
}

// Matrix times vector

std::vector<double> MatrixTimesVector(const std::vector<std::vector<double>>& A, const std::vector<double>& x)
{
	std::vector<double> b(x.size(), 0.0);

	for (size_t i = 0; i < A.size(); i++)
	{
		for (size_t j = 0; j < x.size(); j++)
		{
			b[i] += A[i][j] * x[j];
		}
	}

	return b;
}

// Simple point class

class Point
{
public:

	// Constructors

	Point() : x(0.0), y(0.0)
	{
	}

	Point(double xSource, double ySource) : x(xSource), y(ySource)
	{
	}

	Point(const Point& sourcePoint)
	{
		*this = sourcePoint;
	}

	Point& operator=(const Point& sourcePoint)
	{
		if (this != &sourcePoint)
		{
			x = sourcePoint.x;
			y = sourcePoint.y;
		}

		return *this;
	}

	// Destructor

	~Point()
	{
	}

	// Data members public for convenience

	double x;
	double y;
};

// Simple data point class

class DataPoint
{
public:

	// Constructors

	DataPoint() : x(0.0), y(0.0), value(0.0)
	{
	}

	DataPoint(double xSource, double ySource, double valueSource) : x(xSource), y(ySource), value(valueSource)
	{
	}

	DataPoint(const DataPoint& sourceDataPoint)
	{
		*this = sourceDataPoint;
	}

	DataPoint& operator=(const DataPoint& sourceDataPoint)
	{
		if (this != &sourceDataPoint)
		{
			x = sourceDataPoint.x;
			y = sourceDataPoint.y;
			value = sourceDataPoint.value;
		}

		return *this;
	}

	// Destructor

	~DataPoint()
	{
	}

	// Data members public for convenience

	double x;
	double y;

	double value;
};

// Simple Raster class

class _RasterContext
{
public:

	// Constructor

	_RasterContext::_RasterContext()
	{
	}

	_RasterContext(double xOffsetSource, double yOffsetSource,
		double xScaleSource, double yScaleSource,
		unsigned int widthSource, unsigned int heightSource, bool yFlippedSource = true) : xOffset(xOffsetSource), yOffset(yOffsetSource),
		xScale(xScaleSource), yScale(yScaleSource),
		width(widthSource), height(heightSource), yFlipped(yFlippedSource)
	{
	}

	// Copy constructor

	_RasterContext(const _RasterContext& sourceRasterContext)
	{
		*this = sourceRasterContext;
	}

	// Destructor

	~_RasterContext()
	{
	}

	// Assignment operator

	_RasterContext& operator=(const _RasterContext& sourceRasterContext)
	{
		if (this != &sourceRasterContext)
		{
			xOffset = sourceRasterContext.xOffset;
			yOffset = sourceRasterContext.yOffset;

			xScale = sourceRasterContext.xScale;
			yScale = sourceRasterContext.yScale;

			width = sourceRasterContext.width;
			height = sourceRasterContext.height;

			yFlipped = sourceRasterContext.yFlipped;
		}

		return *this;
	}

	Point XYtoPoint(unsigned int x, unsigned int y)
	{
		double xCoordinate = 0.0;
		double yCoordinate = 0.0;

		xCoordinate = (x * xScale) + xOffset;

		if (yFlipped)
		{
			yCoordinate = ((height - y) * yScale) + yOffset;
		}
		else
		{
			yCoordinate = (y * yScale) + yOffset;
		}

		return Point(xCoordinate, yCoordinate);
	}

	// Data members, public for convenience

	double xOffset = 0.0;
	double yOffset = 0.0;

	double xScale = 0.0;
	double yScale = 0.0;

	unsigned int width = 0;
	unsigned int height = 0;

	bool yFlipped = true;
};

// Cholesky matrix decomposition to lower triangular matrix and its conjugate transpose
// 
// Restricted to positive-definite matrices

class CholeskyDecomposition
{
public:

	// Constructor

	// Matrix is decomposed in-place

	CholeskyDecomposition(std::vector<std::vector<double>>& sourceMatrix);

	// Destructor

	~CholeskyDecomposition();

	// Decomposition into triangular matrices

	bool Decompose();

	// Solve for x in form Ax = b.  A is the original input matrix.

	std::vector<double> Solve(const std::vector<double>& b);

protected:

	CholeskyDecomposition(const CholeskyDecomposition&) = delete;

	void operator=(const CholeskyDecomposition&) = delete;

	// Input matrix

	std::vector<std::vector<double>>& decomposedMatrix;
};

// Constructor

CholeskyDecomposition::CholeskyDecomposition(std::vector<std::vector<double>>& sourceMatrix) : decomposedMatrix(sourceMatrix)
{
	assert(sourceMatrix.size() > 0 && sourceMatrix.size() == sourceMatrix[0].size());
}

// Destructor

CholeskyDecomposition::~CholeskyDecomposition()
{
}

// Decomposition into triangular matrices

bool CholeskyDecomposition::Decompose()
{
	// Enumerate matrix columnwise

	for (size_t j = 0; j < decomposedMatrix.size(); j++)
	{
		for (size_t i = j; i < decomposedMatrix.size(); i++)
		{
			if (i == j)
			{
				double sum = 0.0;

				for (size_t k = 0; k < i; k++)
				{
					sum += std::pow(decomposedMatrix[i][k], 2.0);
				}

				if (decomposedMatrix[i][j] - sum <= 0.0)
				{
					// Not positive definite matrix

					return false;
				}

				decomposedMatrix[i][j] = std::sqrt(decomposedMatrix[i][j] - sum);
			}
			else
			{
				double sum = 0.0;

				for (size_t k = 0; k < j; k++)
				{
					sum += (decomposedMatrix[i][k] * decomposedMatrix[j][k]);
				}

				decomposedMatrix[i][j] = (1 / decomposedMatrix[j][j]) * (decomposedMatrix[i][j] - sum);
				decomposedMatrix[j][i] = decomposedMatrix[i][j];
			}
		}
	}

	return true;
}

// Solve for x in form Ax = b.  A is the original input matrix.

std::vector<double> CholeskyDecomposition::Solve(const std::vector<double>& b)
{
	std::vector<double> y(b.size());

	// First solve lower triangular * y = b with forward substitution

	for (size_t i = 0; i < b.size(); i++)
	{
		double sum = 0.0;

		for (size_t j = 0; j < i; j++)
		{
			sum += (decomposedMatrix[i][j] * y[j]);
		}

		y[i] = (b[i] - sum) / decomposedMatrix[i][i];
	}

	// Now solve upper triangular (transpose of lower triangular) * x = y with back substitution.
	// Note that x can be solved in place using the existing y vector.  No need to allocate 
	// another vector.

	for (int i = static_cast<int>(b.size()) - 1; i >= 0; i--)
	{
		double sum = 0.0;

		for (int j = static_cast<int>(b.size()) - 1; j > i; j--)
		{
			sum += (decomposedMatrix[i][j] * y[j]);
		}

		y[i] = (y[i] - sum) / decomposedMatrix[i][i];
	}

	return y;
}

// Lower/upper decomposition of matrix into a lower triangular matrix and a upper triangular matrix.

class LUDecomposition
{
public:

	// Constructor

	// Matrix is decomposed in-place

	LUDecomposition(std::vector<std::vector<double>>& sourceMatrix);

	// Destructor

	~LUDecomposition();

	// Decomposition into triangular matrices

	bool Decompose();

	// Solve for x in form Ax = b.  A is the original input matrix.

	std::vector<double> Solve(const std::vector<double>& b);

protected:

	LUDecomposition(const LUDecomposition&) = delete;

	void operator=(const LUDecomposition&) = delete;

	// Output matrix after decomposition

	std::vector<std::vector<double>>& decomposedMatrix;

	// Permutation of rows during pivoting

	std::vector<int> rowPermutation;
};

// Constructor

LUDecomposition::LUDecomposition(std::vector<std::vector<double>>& sourceMatrix) : decomposedMatrix(sourceMatrix)
{
	assert(sourceMatrix.size() > 0 && sourceMatrix.size() == sourceMatrix[0].size());
}

// Destructor

LUDecomposition::~LUDecomposition()
{
}

// Decomposition into triangular matrices

bool LUDecomposition::Decompose()
{
	// Initialize the permutation vector

	size_t n = decomposedMatrix.size();

	rowPermutation.reserve(n);

	for (size_t i = 0; i < n; i++)
	{
		rowPermutation.push_back((int)i);
	}

	double det = 1.0;

	// LU factorization.

	for (size_t p = 1; p <= n - 1; p++)
	{
		// Find pivot element.

		for (size_t i = p + 1; i <= n; i++)
		{
			if (std::fabs(decomposedMatrix[rowPermutation[i - 1]][p - 1]) > std::fabs(decomposedMatrix[rowPermutation[p - 1]][p - 1]))
			{
				// Switch the index for the p-1 pivot row if necessary.

				std::swap(rowPermutation[p - 1], rowPermutation[i - 1]);

				det = -det;
			}
		}

		if (decomposedMatrix[rowPermutation[p - 1]][p - 1] == 0.0)
		{
			// The matrix is singular, at least to precision of algorithm

			return false;
		}

		// Multiply the diagonal elements.

		det = det * decomposedMatrix[rowPermutation[p - 1]][p - 1];

		// Form multiplier.

		for (size_t i = p + 1; i <= n; i++)
		{
			decomposedMatrix[rowPermutation[i - 1]][p - 1] /= decomposedMatrix[rowPermutation[p - 1]][p - 1];

			// Eliminate [p-1].

			for (size_t j = p + 1; j <= n; j++)
			{
				decomposedMatrix[rowPermutation[i - 1]][j - 1] -= decomposedMatrix[rowPermutation[i - 1]][p - 1] * decomposedMatrix[rowPermutation[p - 1]][j - 1];
			}
		}
	}

	det = det * decomposedMatrix[rowPermutation[n - 1]][n - 1];

	return (det != 0.0);
}

// Solve for x in form Ax = b.  A is the original input matrix.

// Note: b is modified in-place for row permutations

std::vector<double> LUDecomposition::Solve(const std::vector<double>& b)
{
	// Our decomposed matrix is comprised of both the lower and upper diagonal matrices.

	// The rows of this matrix have been permutated during the decomposition process.  The
	// rowPermutation indicates the proper row order.

	// The lower diagonal matrix only include elements below the diagonal with diagonal 
	// elements set to 1.

	// The upper diagonal matrix is fully specified.

	// First solve Ly = Pb for y using forward substitution. P is a permutated identity matrix.

	std::vector<double> y(b.size());

	for (size_t i = 0; i < y.size(); i++)
	{
		size_t currentRow = rowPermutation[i];

		double sum = 0.0;

		for (size_t j = 0; j < i; j++)
		{
			sum += (decomposedMatrix[currentRow][j] * y[j]);
		}

		y[i] = (b[currentRow] - sum);
	}

	// Now solve Ux = y for x using back substitution.  Note that 
	// x can be solved in place using the existing y vector.  No need
	// to allocate another vector.

	for (int i = static_cast<int>(b.size()) - 1; i >= 0; i--)
	{
		size_t currentRow = rowPermutation[i];

		double sum = 0.0;

		for (int j = static_cast<int>(b.size()) - 1; j > i; j--)
		{
			sum += (decomposedMatrix[currentRow][j] * y[j]);
		}

		y[i] = (y[i] - sum) / decomposedMatrix[currentRow][i];
	}

	return y;
}

// Simple and Ordinary Kriging

class _Kriging
{
public:

	// Data model

	enum _Model
	{
		Linear,
		LinearWithoutIntercept,
		Spherical,
		Exponential,
		Gaussian
	};

	// Constructor

	_Kriging(std::vector<DataPoint>& dataPoint, _RasterContext& rasterContext);

	// Destructor

	~_Kriging();

	// Load the input data

	void Initialize();

	// Calculate the semivariances

	void CalculateExperimentalVariogram(double lag, double lagTolerance);

	// Calculate models

	std::vector<double> CalculateLinearModel() const;
	std::vector<double> CalculateLinearModelWithoutIntercept(double nugget) const;
	std::vector<double> CalculateSphericalModel(double nugget, double sill, double range) const;
	std::vector<double> CalculateExponentialModel(double nugget, double sill, double range) const;
	std::vector<double> CalculateGaussianModel(double nugget, double sill, double range) const;

	// Retrieve semivariogram parameters

	double GetEstimatedNugget() const
	{
		return estimatedNugget;
	}

	double GetEstimatedRange() const
	{
		return estimatedRange;
	}

	double GetEstimatedSill() const
	{
		return estimatedSill;
	}

	// Get lag parameters

	double GetDefaultLag() const
	{
		return lag;
	}

	double GetDefaultLagTolerance() const
	{
		return lagTolerance;
	}

	// Get empirical semivariogram vectors

	const std::vector<double>& GetLagDistances()
	{
		return lagDistance;
	}

	const std::vector<double>& GetLagSemivariances()
	{
		return lagSemivariance;
	}

	const std::vector<unsigned int>& GetLagCounts()
	{
		return lagCount;
	}

	// Simple Kriging

	void SimpleKrige(_Model model, double nugget, double sill, double range);
	void SimpleKrige(_Model model,void* buffer, double nugget, double sill, double range);

	// Ordinary Kriging

	void OrdinaryKrige(_Model model, double nugget, double sill, double range,
		unsigned int minPoints, unsigned int maxPoints, double maxDistance);
	void OrdinaryKrige(_Model model, void* buffer, double nugget, double sill, double range,
		unsigned int minPoints, unsigned int maxPoints, double maxDistance);

public:

	// Fill a map of all distances with semivariogram

	void CalculateDistanceMap();

	// Fill a matrix of variograms over all point distances

	std::vector<std::vector<double>> CalculateVariogramMatrix(const std::vector<DataPoint>& dataPointCandidate,
		_Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const;

	// Fill a matrix of covariograms over all point distances

	std::vector<std::vector<double>> CalculateCovariogramMatrix(const std::vector<DataPoint>& dataPointCandidate,
		_Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const;

	// Fill a vector of variograms over all distances from a given point

	std::vector<double> CalculateVariogramVector(const std::vector<DataPoint>& dataPointCandidate,
		double xCoordinate, double yCoordinate, _Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const;

	// Fill a vector of covariograms over all distances from a given point

	std::vector<double> CalculateCovariogramVector(const std::vector<DataPoint>& dataPointCandidate,
		double xCoordinate, double yCoordinate, _Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const;

	// Fill a map of distances from a fixed point with model variogram

	std::multimap<double, size_t> CalculateDistanceMapForPoint(double pointX, double pointY, unsigned int maxPoints, double maxDistance) const;

	// Find the appropriate lag parameters

	void CalculateDefaultLagParameters();

	// Find the estimated sill (sample variance)

	void CalculateEstimatedVariogramParameters();

	// Calculate the variogram

	double CalculateVariogram(_Model model, double distance, double nugget, double sill, double range) const;

	// Calculate the covariogram

	double CalculateCovariogram(_Model model, double distance, double nugget, double sill, double range) const;

	// Simple linear regression

	void DoSimpleLinearRegression(const std::vector<double>& X, const std::vector<double>& Y, double* slope, double* intercept) const;
	void DoSimpleLinearRegressionWithoutIntercept(const std::vector<double>& X, const std::vector<double>& Y, double* slope, double intercept) const;

	// Simple Kringing for an individual point

	double SimpleKrigeForPoint(double xCoordinate, double yCoordinate,
		_Model model, double nugget, double sill, double range,
		CholeskyDecomposition& cholesky, std::vector<double>& residuals, double estimatedMean);

	// Ordinary Kringing for an individual point

	double OrdinaryKrigeForPoint(double xCoordinate, double yCoordinate,
		_Model model, double nugget, double sill, double range,
		LUDecomposition& luDecomposition, const std::vector<DataPoint>& dataPointCandidate);

	double estimatedNugget = 0.0;
	double estimatedSill = 0.0;
	double estimatedRange = 0.0;

	double lag = 0.0;
	double lagTolerance = 0.0;

	std::vector<DataPoint> dataPoint;

	std::multimap<double, double> semiVariogram;

	std::vector<double> lagDistance;
	std::vector<double> lagSemivariance;
	std::vector<unsigned int> lagCount;

	_RasterContext rasterContext;
};

// Constructor

_Kriging::_Kriging(std::vector<DataPoint>& sourceDataPoint, _RasterContext& rasterContext) : rasterContext(rasterContext)
{
	dataPoint.swap(sourceDataPoint);
}

// Intialize internals

void _Kriging::Initialize()
{
	CalculateDistanceMap();
	CalculateDefaultLagParameters();
	CalculateExperimentalVariogram(GetDefaultLag(), GetDefaultLagTolerance());
	CalculateEstimatedVariogramParameters();
}

// Destructor

_Kriging::~_Kriging()
{
}

// Simple Kriging

void _Kriging::SimpleKrige(_Model model, double nugget, double sill, double range)
{
#if 0

	// Simplistic test data for one point

	sill = 0.77757733618823532;

	// x = { 2700, 2300, 900, 900, 500, 3700 };
	// y = { 4300, 5700, 5100, 3700, 4900, 5100 };
	// z = { 12.149, 12.681, 14.414, 13.835, 14.591, 12.867 };

	dataPoint.resize(6);

	dataPoint[0] = DataPoint(2700, 4300, 12.149);
	dataPoint[1] = DataPoint(2300, 5700, 12.681);
	dataPoint[2] = DataPoint(900, 5100, 14.414);
	dataPoint[3] = DataPoint(900, 3700, 13.835);
	dataPoint[4] = DataPoint(500, 4900, 14.591);
	dataPoint[5] = DataPoint(3700, 5100, 12.867);

	// Find distances between all points and calculate covariograms

	std::vector<std::vector<double>> distanceCovariogramMatrix = CalculateCovariogramMatrix(dataPoint, model, nugget, sill, range, false);

	// Decompose the covariogram matrix 

	CholeskyDecomposition choleskyDecomposition(distanceCovariogramMatrix);

	choleskyDecomposition.Decompose();

	// Calculate mean and residuals.

	// Note: For this example, we are using the mean over the entire sample data set and not the mean of the datapoints.

	double estimatedMean = 14.695879999999999;

	std::vector<double> residuals(dataPoint.size());

	std::transform(dataPoint.begin(), dataPoint.end(), residuals.begin(), [&](const DataPoint& current) { return current.value - estimatedMean; });

	// Should be approximately 12.826384946039500...

	double estimatedPoint = SimpleKrigeForPoint(2000, 4700, model, nugget, sill, range, choleskyDecomposition, residuals, estimatedMean);

	return;

#endif

#if 1

	// Find distances between all points and calculate covariograms

	std::vector<std::vector<double>> distanceCovariogramMatrix = CalculateCovariogramMatrix(dataPoint, model, nugget, sill, range, false);

	// Decompose the covariogram matrix 

	CholeskyDecomposition choleskyDecomposition(distanceCovariogramMatrix);

	choleskyDecomposition.Decompose();

	// Estimate mean and calculate residuals

	double estimatedMean = std::accumulate(dataPoint.begin(), dataPoint.end(), 0.0, [](double sum, const DataPoint& dataPoint) { return sum + dataPoint.value; }) / dataPoint.size();

	std::vector<double> residuals(dataPoint.size());

	std::transform(dataPoint.begin(), dataPoint.end(), residuals.begin(), [&](const DataPoint& current) { return current.value - estimatedMean; });

	double estimatedZ = 0.0;

	// For now, output to standard out to plot in Excel.  Later we will just assign the kriged estimate to the data grid.

	// Y scale across row 1

	for (unsigned int j = 0; j < rasterContext.height; j++)
	{
		std::cout << ", " << rasterContext.XYtoPoint(0, j).y;
	}

	std::cout << std::endl;

	for (unsigned int i = 0; i < rasterContext.width; i++)
	{
		// X scale down column 1

		std::cout << rasterContext.XYtoPoint(i, 0).x;

		for (unsigned int j = 0; j < rasterContext.height; j++)
		{
			// The kriged estimate
			//估计
			Point point = rasterContext.XYtoPoint(i, j);

			estimatedZ = SimpleKrigeForPoint(point.x, point.y, model, nugget, sill, range, choleskyDecomposition, residuals, estimatedMean);

			std::cout << ", " << estimatedZ;
		}

		std::cout << std::endl;
	}

	return;

#endif
}

void _Kriging::SimpleKrige(_Model model, void* buf, double nugget, double sill, double range)
{
	float* buffer = static_cast<float*>(buf);
	// Find distances between all points and calculate covariograms

	std::vector<std::vector<double>> distanceCovariogramMatrix = CalculateCovariogramMatrix(dataPoint, model, nugget, sill, range, false);

	// Decompose the covariogram matrix 

	CholeskyDecomposition choleskyDecomposition(distanceCovariogramMatrix);

	choleskyDecomposition.Decompose();

	// Estimate mean and calculate residuals

	double estimatedMean = std::accumulate(dataPoint.begin(), dataPoint.end(), 0.0, [](double sum, const DataPoint& dataPoint) { return sum + dataPoint.value; }) / dataPoint.size();

	std::vector<double> residuals(dataPoint.size());

	std::transform(dataPoint.begin(), dataPoint.end(), residuals.begin(), [&](const DataPoint& current) { return current.value - estimatedMean; });

	double estimatedZ = 0.0;

	// For now, output to standard out to plot in Excel.  Later we will just assign the kriged estimate to the data grid.

	// Y scale across row 1

	//for (unsigned int j = 0; j < rasterContext.height; j++)
	//{
	//	std::cout << ", " << rasterContext.XYtoPoint(0, j).y;
	//}

	//std::cout << std::endl;

	for (unsigned int i = 0; i < rasterContext.width; i++)
	{
		// X scale down column 1

		//std::cout << rasterContext.XYtoPoint(i, 0).x;

		for (unsigned int j = 0; j < rasterContext.height; j++)
		{
			// The kriged estimate
			//估计
			Point point = rasterContext.XYtoPoint(i, j);

			estimatedZ = SimpleKrigeForPoint(point.x, point.y, model, nugget, sill, range, choleskyDecomposition, residuals, estimatedMean);
			buffer[j * rasterContext.width + i] = estimatedZ;
			//std::cout << ", " << estimatedZ;
		}

		//std::cout << std::endl;
	}

	return;
}

// Ordinary Kriging

void _Kriging::OrdinaryKrige(_Model model, double nugget, double sill, double range,
	unsigned int minPoints, unsigned int maxPoints, double maxDistance)
{
#if 0

	// Simplistic test data for one point

	sill = 0.77757733618823532;

	// x = { 2700, 2300, 900, 900, 500, 3700 };
	// y = { 4300, 5700, 5100, 3700, 4900, 5100 };
	// z = { 12.149, 12.681, 14.414, 13.835, 14.591, 12.867 };

	dataPoint.resize(6);

	dataPoint[0] = DataPoint(2700, 4300, 12.149);
	dataPoint[1] = DataPoint(2300, 5700, 12.681);
	dataPoint[2] = DataPoint(900, 5100, 14.414);
	dataPoint[3] = DataPoint(900, 3700, 13.835);
	dataPoint[4] = DataPoint(500, 4900, 14.591);
	dataPoint[5] = DataPoint(3700, 5100, 12.867);

	// Find distances between all points and calculate covariograms

	std::vector<std::vector<double>> distanceCovariogramMatrix = CalculateCovariogramMatrix(dataPoint, model, nugget, sill, range, true);

	// Decompose the covariogram matrix.  Because of the Lagrange multiplier, this will not be a positive definite matrix and we will
	// need to use LU decomposition.

	LUDecomposition luDecomposition(distanceCovariogramMatrix);

	luDecomposition.Decompose();

	// Should be approximately 12.931105462857452...

	double estimatedPoint = OrdinaryKrigeForPoint(2000, 4700, model, nugget, sill, range, luDecomposition, dataPoint);

	return;

#endif

#if 1

	// For now, output to standard out to plot in Excel.  Later will be just assign the kriged estimate to the data grid.

	// Y scale across row 1

	for (unsigned int j = 0; j < rasterContext.height; j++)
	{
		std::cout << ", " << rasterContext.XYtoPoint(0, j).y;
	}

	std::cout << std::endl;

	std::vector<DataPoint> dataPointCandidate(maxPoints);

	if (maxDistance == 0.0)
	{
		maxDistance = std::numeric_limits<double>().max();
	}

	for (unsigned int i = 0; i < rasterContext.width; i++)
	{
		// X scale down column 1

		std::cout << rasterContext.XYtoPoint(i, 0).x;

		for (unsigned int j = 0; j < rasterContext.height; j++)
		{
			// Current x, y values

			Point point = rasterContext.XYtoPoint(i, j);

			// Find our candidate points.

			std::multimap<double, size_t> candidateDistanceMap = CalculateDistanceMapForPoint(point.x, point.y, maxPoints, maxDistance);

			double estimatedZ = 0.0;

			if (candidateDistanceMap.size() == minPoints)
			{
				dataPointCandidate.resize(candidateDistanceMap.size());

				size_t index = 0;

				for (std::multimap<double, size_t>::const_iterator iterator = candidateDistanceMap.begin();
					iterator != candidateDistanceMap.end();
					iterator++)
				{
					dataPointCandidate[index] = this->dataPoint[iterator->second];

					index++;
				}

				// Find distances between all points and calculate covariograms

				std::vector<std::vector<double>> distanceCovariogramMatrix = CalculateCovariogramMatrix(dataPointCandidate, model, nugget, sill, range, true);

				// Decompose the covariogram matrix.  Because of the Lagrange multiplier, this will not be a positive definite matrix and we will
				// need to use LU decomposition.

				LUDecomposition luDecomposition(distanceCovariogramMatrix);

				luDecomposition.Decompose();

				estimatedZ = OrdinaryKrigeForPoint(point.x, point.y, model, nugget, sill, range, luDecomposition, dataPointCandidate);
			}

			std::cout << ", " << estimatedZ;
		}

		std::cout << std::endl;
	}

	return;

#endif
}

void _Kriging::OrdinaryKrige(_Model model, void* buf, double nugget, double sill, double range, unsigned int minPoints, unsigned int maxPoints, double maxDistance)
{
	float* buffer = static_cast<float*>(buf);
	// For now, output to standard out to plot in Excel.  Later will be just assign the kriged estimate to the data grid.

	// Y scale across row 1

	//for (unsigned int j = 0; j < rasterContext.height; j++)
	//{
	//	std::cout << ", " << rasterContext.XYtoPoint(0, j).y;
	//}

	//std::cout << std::endl;

	std::vector<DataPoint> dataPointCandidate(maxPoints);

	if (maxDistance == 0.0)
	{
		maxDistance = std::numeric_limits<double>().max();
	}

	for (unsigned int i = 0; i < rasterContext.width; i++)
	{
		// X scale down column 1

		//std::cout << rasterContext.XYtoPoint(i, 0).x;

		for (unsigned int j = 0; j < rasterContext.height; j++)
		{
			// Current x, y values

			Point point = rasterContext.XYtoPoint(i, j);

			// Find our candidate points.

			std::multimap<double, size_t> candidateDistanceMap = CalculateDistanceMapForPoint(point.x, point.y, maxPoints, maxDistance);

			double estimatedZ = 0.0;

			if (candidateDistanceMap.size() >= minPoints)
			{
				dataPointCandidate.resize(candidateDistanceMap.size());

				size_t index = 0;

				for (std::multimap<double, size_t>::const_iterator iterator = candidateDistanceMap.begin();
					iterator != candidateDistanceMap.end();
					iterator++)
				{
					dataPointCandidate[index] = this->dataPoint[iterator->second];

					index++;
				}

				// Find distances between all points and calculate covariograms

				std::vector<std::vector<double>> distanceCovariogramMatrix = CalculateCovariogramMatrix(dataPointCandidate, model, nugget, sill, range, true);

				// Decompose the covariogram matrix.  Because of the Lagrange multiplier, this will not be a positive definite matrix and we will
				// need to use LU decomposition.

				LUDecomposition luDecomposition(distanceCovariogramMatrix);

				luDecomposition.Decompose();

				estimatedZ = OrdinaryKrigeForPoint(point.x, point.y, model, nugget, sill, range, luDecomposition, dataPointCandidate);
			}
			buffer[j * rasterContext.width + i] = estimatedZ;
			//std::cout << ", " << estimatedZ;
		}

		//std::cout << std::endl;
	}

	return;
}

// Simple Kringing for an individual point

double _Kriging::SimpleKrigeForPoint(double xCoordinate, double yCoordinate,
	_Model model, double nugget, double sill, double range,
	CholeskyDecomposition& choleskyDecomposition, std::vector<double>& residuals, double estimatedMean)
{
	// Find distances over given points and calculate covariograms

	std::vector<double> distanceCovariogramVector = CalculateCovariogramVector(dataPoint, xCoordinate, yCoordinate, model, nugget, sill, range, false);

	// Solve Ax = b, for x.  A represents the covariogram matrix of all distances and b is the covariogram vector for the current point.
	// x will be a vector of weights.

	std::vector<double> weights = choleskyDecomposition.Solve(distanceCovariogramVector);

	// Multiply the weights by the residuals and add the estimated mean to yield estimate for this point.

	return std::inner_product(weights.begin(), weights.end(), residuals.begin(), 0.0) + estimatedMean;
}

// Ordinary Kringing for an individual point

double _Kriging::OrdinaryKrigeForPoint(double xCoordinate, double yCoordinate,
	_Model model, double nugget, double sill, double range,
	LUDecomposition& luDecomposition, const std::vector<DataPoint>& dataPointCandidate)
{
	// Find distances over given points and calculate covariograms

	std::vector<double> distanceCovariogramVector = CalculateCovariogramVector(dataPointCandidate, xCoordinate, yCoordinate, model, nugget, sill, range, true);

	// Solve Ax = b, for x.  A represents the covariogram matrix of all distances and b is the covariogram vector for the current point.
	// x will be a vector of weights.

	std::vector<double> weights = luDecomposition.Solve(distanceCovariogramVector);

	// Multiply the weights by the residuals to yield estimate for this point.

	return std::inner_product(weights.begin(), weights.end() - 1, dataPointCandidate.begin(), 0.0, std::plus<double>(), [](double weight, const DataPoint& dataPointCandidate) { return weight * dataPointCandidate.value; });
}

// Fill a map of distances from a fixed point with variogram

std::multimap<double, size_t> _Kriging::CalculateDistanceMapForPoint(double pointX, double pointY, unsigned int maxPoints, double maxDistance) const
{
	std::multimap<double, size_t> pointDistanceMap;

	double currentDistance = 0.0;

	for (size_t i = 0; i < dataPoint.size(); i++)
	{
		currentDistance = std::sqrt(std::pow(dataPoint[i].x - pointX, 2.0) + std::pow(dataPoint[i].y - pointY, 2.0));

		if (currentDistance <= maxDistance)
		{
			if (pointDistanceMap.size() == maxPoints)
			{
				if (currentDistance < pointDistanceMap.rbegin()->first)
				{
					pointDistanceMap.insert(std::pair<double, size_t>(currentDistance, i));

					pointDistanceMap.erase(pointDistanceMap.rbegin()->first);
				}
			}
			else
			{
				pointDistanceMap.insert(std::pair<double, size_t>(currentDistance, i));
			}
		}
	}

	return pointDistanceMap;
}

// Fill a vector of variograms over all distances from a given point

std::vector<double> _Kriging::CalculateVariogramVector(const std::vector<DataPoint>& dataPointCandidate,
	double xCoordinate, double yCoordinate, _Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const
{
	std::vector<double> distanceVector(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), LagrangeMultiplier ? 1.0 : 0.0);

	double distance = 0.0;
	double variogram = 0.0;

	for (size_t i = 0; i < dataPointCandidate.size(); i++)
	{
		distance = std::sqrt(std::pow(dataPointCandidate[i].x - xCoordinate, 2.0) + std::pow(dataPointCandidate[i].y - yCoordinate, 2.0));

		variogram = CalculateVariogram(model, distance, nugget, sill, range);

		distanceVector[i] = variogram;
	}

	return distanceVector;
}

// Fill a vector of covariograms over all distances from a given point

std::vector<double> _Kriging::CalculateCovariogramVector(const std::vector<DataPoint>& dataPointCandidate,
	double xCoordinate, double yCoordinate, _Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const
{
	std::vector<double> distanceVector(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), LagrangeMultiplier ? 1.0 : 0.0);

	double distance = 0.0;
	double covariogram = 0.0;

	for (size_t i = 0; i < dataPointCandidate.size(); i++)
	{
		distance = std::sqrt(std::pow(dataPointCandidate[i].x - xCoordinate, 2.0) + std::pow(dataPointCandidate[i].y - yCoordinate, 2.0));

		covariogram = CalculateCovariogram(model, distance, nugget, sill, range);

		distanceVector[i] = covariogram;
	}

	return distanceVector;
}

// Fill a matrix of variograms over all point distances

std::vector<std::vector<double>> _Kriging::CalculateVariogramMatrix(const std::vector<DataPoint>& dataPointCandidate,
	_Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const
{
	std::vector<std::vector<double>> distanceMatrix(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), std::vector<double>(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), LagrangeMultiplier ? 1.0 : 0.0));

	double distance = 0.0;
	double variogram = 0.0;

	for (size_t i = 0; i < dataPointCandidate.size(); i++)
	{
		distanceMatrix[i][i] = CalculateVariogram(model, 0.0, nugget, sill, range);

		for (size_t j = i + 1; j < dataPointCandidate.size(); j++)
		{
			distance = std::sqrt(std::pow(dataPointCandidate[i].x - dataPointCandidate[j].x, 2.0) + std::pow(dataPointCandidate[i].y - dataPointCandidate[j].y, 2.0));

			variogram = CalculateVariogram(model, distance, nugget, sill, range);

			distanceMatrix[i][j] = variogram;
			distanceMatrix[j][i] = distanceMatrix[i][j];
		}
	}

	if (LagrangeMultiplier)
	{
		distanceMatrix[dataPointCandidate.size()][dataPointCandidate.size()] = 0.0;
	}

	return distanceMatrix;
}

// Fill a matrix of covariograms over all point distances

std::vector<std::vector<double>> _Kriging::CalculateCovariogramMatrix(const std::vector<DataPoint>& dataPointCandidate,
	_Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const
{
	std::vector<std::vector<double>> distanceMatrix(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), std::vector<double>(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), LagrangeMultiplier ? 1.0 : 0.0));

	double distance = 0.0;
	double covariogram = 0.0;

	for (size_t i = 0; i < dataPointCandidate.size(); i++)
	{
		distanceMatrix[i][i] = CalculateCovariogram(model, 0.0, nugget, sill, range);

		for (size_t j = i + 1; j < dataPointCandidate.size(); j++)
		{
			distance = std::sqrt(std::pow(dataPointCandidate[i].x - dataPointCandidate[j].x, 2.0) + std::pow(dataPointCandidate[i].y - dataPointCandidate[j].y, 2.0));

			covariogram = CalculateCovariogram(model, distance, nugget, sill, range);

			distanceMatrix[i][j] = covariogram;
			distanceMatrix[j][i] = distanceMatrix[i][j];
		}
	}

	if (LagrangeMultiplier)
	{
		distanceMatrix[dataPointCandidate.size()][dataPointCandidate.size()] = 0.0;
	}

	return distanceMatrix;
}

// Fill a map of all distances with semivariogram

void _Kriging::CalculateDistanceMap()
{
	// Rather than take the simple variance of the variable, a better estimation should be the
	// average variance amoung all distances.

	estimatedSill = 0.0;

	double variance = 0.0;

	for (size_t i = 0; i < dataPoint.size(); i++)
	{
		for (size_t j = i + 1; j < dataPoint.size(); j++)
		{
			variance = std::pow(dataPoint[i].value - dataPoint[j].value, 2.0);

			estimatedSill += variance;

			semiVariogram.insert(std::pair<double, double>(std::sqrt(std::pow(dataPoint[i].x - dataPoint[j].x, 2.0) + std::pow(dataPoint[i].y - dataPoint[j].y, 2.0)), variance));
		}
	}

	estimatedSill = (estimatedSill / 2.0) / ((dataPoint.size() * (dataPoint.size() - 1)) / 2.0);
}

// Find the appropriate lag parameters

void _Kriging::CalculateDefaultLagParameters()
{
	// This algorithm seems to come up with reasonable lag parameters.

	double minDistance = semiVariogram.begin()->first;

	lag = minDistance * 2 * 1.25;

	double lagBoundary = (semiVariogram.rbegin()->first - semiVariogram.begin()->first) / 2.0;

	unsigned int lagCount = static_cast<unsigned int>(lagBoundary / lag + 0.5);

	const unsigned MINIMUM_LAG_COUNT = 20;

	if (lagCount < MINIMUM_LAG_COUNT)
	{
		lag = minDistance;
	}

	lagTolerance = static_cast<unsigned int>((lag / 2.0) + 0.5);
}

// Find the estimated variogram parameters

void _Kriging::CalculateEstimatedVariogramParameters()
{
	// Note: Determination of sill and range are only rough estimates.
	//       It is up to the user to interpret the empirical and model
	//       variogram plots to assess the validity of the parameters 
	//       used in the chosen model.

	// Rather than take the simple variance of the variable, a better estimation should be the
	// average variance amoung all distances.

	// Sill is the simple variance

	//double mean = std::accumulate(z.begin(), z.end(), 0.0) / z.size();
	//double variance = 0.0;

	//std::for_each(std::begin(z), std::end(z), [&](const double d)
	//{
	//	variance += (d - mean) * (d - mean);
	//});

	//sill = variance / (z.size() - 1);

	// For fixed models, range is the first distance where the sill is reached.

	// For asymptotic models it would be the first distance where the semivariance
	// reaches 95% of the sill.

	// We will assume a fixed model to start.

	estimatedRange = 0.0;

	std::vector<double>::const_iterator semivariance;
	std::vector<double>::const_iterator distance;

	for (semivariance = lagSemivariance.begin(), distance = lagDistance.begin(); semivariance != lagSemivariance.end(); semivariance++, distance++)
	{
		if (*semivariance > estimatedSill)
		{
			estimatedRange = (estimatedRange + *distance) / 2.0;
			break;
		}
		else
		{
			estimatedRange = *distance;
		}
	}

	// Calculate estimated nugget

	double slope = 0.0;	// not used

	DoSimpleLinearRegression(lagDistance, lagSemivariance, &slope, &estimatedNugget);
}

// Simple linear regression

void _Kriging::DoSimpleLinearRegression(const std::vector<double>& X, const std::vector<double>& Y, double* slope, double* intercept) const
{
	double Xmean = std::accumulate(X.begin(), X.end(), 0.0) / X.size();
	double Ymean = std::accumulate(Y.begin(), Y.end(), 0.0) / Y.size();

	std::vector<double>::const_iterator Xi = X.begin();
	std::vector<double>::const_iterator Yi = Y.begin();

	double numerator = 0.0;
	double denominator = 0.0;

	while (Xi != lagDistance.end())
	{
		numerator += ((*Xi - Xmean) * (*Yi - Ymean));
		denominator += ((*Xi - Xmean) * (*Xi - Xmean));

		++Xi;
		++Yi;
	}

	*slope = numerator / denominator;
	*intercept = Ymean - (*slope * Xmean);
}

void _Kriging::DoSimpleLinearRegressionWithoutIntercept(const std::vector<double>& X, const std::vector<double>& Y, double* slope, double intercept) const
{
	double Xmean = std::accumulate(X.begin(), X.end(), 0.0) / X.size();
	double Ymean = std::accumulate(Y.begin(), Y.end(), 0.0) / Y.size();

	std::vector<double>::const_iterator Xi = X.begin();
	std::vector<double>::const_iterator Yi = Y.begin();

	double numerator = 0.0;
	double denominator = 0.0;

	while (Xi != lagDistance.end())
	{
		numerator += (*Xi * (*Yi - intercept));
		denominator += (*Xi * *Xi);

		++Xi;
		++Yi;
	}

	*slope = numerator / denominator;
}

// Calculate the semivariances

void _Kriging::CalculateExperimentalVariogram(double lag, double lagTolerance)
{
	// Clear containers from any previous calculations

	lagDistance.clear();
	lagSemivariance.clear();
	lagCount.clear();

	// Only consider points over half the distance.

	double lagBoundary = (semiVariogram.rbegin()->first - semiVariogram.begin()->first) / 2.0;

	double currentLagDistance = lag / 2.0;

	double currentlagSemivariogram = 0.0;

	unsigned int currentLagCount = 0;

	for (std::multimap<double, double>::const_iterator iterator = semiVariogram.begin(); iterator != semiVariogram.end() && currentLagDistance <= lagBoundary; iterator++)
	{
		if (iterator->first > currentLagDistance + lagTolerance)
		{
			if (currentLagCount > 0)
			{
				lagDistance.push_back(currentLagDistance);
				lagSemivariance.push_back((currentlagSemivariogram / currentLagCount) / 2.0);
				lagCount.push_back(currentLagCount);

				currentlagSemivariogram = 0.0;
				currentLagCount = 0;
			}

			currentLagDistance += lag;
		}

		if (iterator->first <= currentLagDistance + lagTolerance)
		{
			currentlagSemivariogram += iterator->second;
			currentLagCount++;
		}
	}
}

// Calculate the variogram

double _Kriging::CalculateVariogram(_Model model, double distance, double nugget, double sill, double range) const
{
	// Notes
	//
	// Linear models do not use nugget, sill or range terminology. For convenience the intercept is passed
	// as the nugget and the intercept as the sill.

	double variogram = 0.0;

	switch (model)
	{
	case Linear:
		variogram = nugget + (sill * distance);
		break;
	case LinearWithoutIntercept:
		variogram = sill * distance;
		break;
	case Spherical:
		if (distance == 0.0)
		{
			variogram = 0.0;
		}
		else if (distance <= range)
		{

			variogram = nugget + (sill - nugget) * ((1.5 * (distance / range) - 0.5 * (std::pow(distance / range, 3.0))));
		}
		else
		{
			variogram = sill;
		}
		break;
	case Exponential:
		if (distance > 0.0)
		{
			variogram = nugget + (sill - nugget) * (1 - std::exp(-(distance / range)));
		}
		else
		{
			variogram = 0.0;
		}
		break;
	case Gaussian:
		if (distance > 0.0)
		{
			variogram = 0.0;
		}
		else
		{
			variogram = nugget + (sill - nugget) * (1 - std::exp(-std::pow(distance / range, 2.0)));
		}
		break;
	default:
		assert(false);
		break;
	}

	return variogram;
}

// Calculate the covariogram

double _Kriging::CalculateCovariogram(_Model model, double distance, double nugget, double sill, double range) const
{
	// Notes
	//
	// As linear models do not have a sill, it is not possible to calculate a covariogram.

	double covariogram = 0.0;

	switch (model)
	{
	case Linear:
	case LinearWithoutIntercept:
		assert(false);
		break;
	case Spherical:
		if (distance == 0.0)
		{
			covariogram = sill;
		}
		else if (distance <= range)
		{
			covariogram = sill * (1 - (1.5 * (distance / range) - 0.5 * (std::pow(distance / range, 3.0))));
		}
		else
		{
			covariogram = 0.0;
		}
		break;
	case Exponential:
		if (distance == 0.0)
		{
			covariogram = sill;
		}
		else
		{
			covariogram = (sill - nugget) * (std::exp(-distance / range));
		}
		break;
	case Gaussian:
		if (distance == 0.0)
		{
			covariogram = sill;
		}
		else
		{
			covariogram = (sill - nugget) * (std::exp(-std::pow(distance / range, 2.0)));
		}
		break;
	default:
		assert(false);
		break;
	}

	return covariogram;
}

// Calculate models

std::vector<double> _Kriging::CalculateLinearModel() const
{
	std::vector<double> variogram;

	variogram.reserve(lagDistance.size());

	// Get slope and intercept

	double slope = 0.0;
	double intercept = 0.0;

	DoSimpleLinearRegression(lagDistance, lagSemivariance, &slope, &intercept);

	for (std::vector<double>::const_iterator distance = lagDistance.begin(); distance < lagDistance.end(); distance++)
	{
		variogram.push_back(CalculateVariogram(Linear, *distance, intercept, slope, 0.0));
	}

	return variogram;
}

std::vector<double> _Kriging::CalculateLinearModelWithoutIntercept(double nugget) const
{
	std::vector<double> variogram;

	variogram.reserve(lagDistance.size());

	// Get slope

	double slope = 0.0;

	DoSimpleLinearRegressionWithoutIntercept(lagDistance, lagSemivariance, &slope, nugget);

	for (std::vector<double>::const_iterator distance = lagDistance.begin(); distance < lagDistance.end(); distance++)
	{
		variogram.push_back(CalculateVariogram(LinearWithoutIntercept, *distance, 0.0, slope, 0.0));
	}

	return variogram;
}

std::vector<double> _Kriging::CalculateSphericalModel(double nugget, double sill, double range) const
{
	std::vector<double> variogram;

	variogram.reserve(lagDistance.size());

	for (std::vector<double>::const_iterator distance = lagDistance.begin(); distance < lagDistance.end(); distance++)
	{
		variogram.push_back(CalculateVariogram(Spherical, *distance, nugget, sill, range));
	}

	return variogram;
}

std::vector<double> _Kriging::CalculateExponentialModel(double nugget, double sill, double range) const
{
	std::vector<double> variogram;

	variogram.reserve(lagDistance.size());

	for (std::vector<double>::const_iterator distance = lagDistance.begin(); distance < lagDistance.end(); distance++)
	{
		variogram.push_back(CalculateVariogram(Exponential, *distance, nugget, sill, range));
	}

	return variogram;
}

std::vector<double> _Kriging::CalculateGaussianModel(double nugget, double sill, double range) const
{
	std::vector<double> variogram;

	variogram.reserve(lagDistance.size());

	for (std::vector<double>::const_iterator distance = lagDistance.begin(); distance < lagDistance.end(); distance++)
	{
		variogram.push_back(CalculateVariogram(Gaussian, *distance, nugget, sill, range));
	}

	return variogram;
}

// Load the input data

bool KrigeZoneAData()
{
	bool result = false;

	std::vector<DataPoint> dataPoint;

	std::ifstream dataFile(".\\data\\ZoneA.csv", std::ifstream::in);

	if (dataFile.good())
	{
		result = true;

		std::string line;

		// Skip header row

		std::getline(dataFile, line);

		while (!dataFile.eof())
		{
			std::getline(dataFile, line);

			std::stringstream lineStream(line);

			std::string cell;

			unsigned int column = 0;

			double x = 0.0;
			double y = 0.0;

			double value = 0.0;

			while (std::getline(lineStream, cell, ',') && column < 4)
			{
				if (column == 0)
				{
					x = std::atof(cell.c_str());
				}
				else if (column == 1)
				{
					y = std::atof(cell.c_str());
				}
				else if (column == 3)
				{
					value = std::atof(cell.c_str());

					dataPoint.push_back(DataPoint(x, y, value));
				}

				column++;
			}
		}
	}

	if (result)
	{
		_RasterContext rasterContext(0, 0, 100, 100, 200, 160, false);

		_Kriging kriging(dataPoint, rasterContext);

		kriging.Initialize();

		//std::cout << "Distance" << ", " << "Semivariance" << ", " << "Counts" << ", " << "Linear Variogram" << std::endl;

		//std::vector<double>::const_iterator distance = kriging.GetLagDistances().begin();
		//std::vector<double>::const_iterator semivariance = kriging.GetLagSemivariances().begin();
		//std::vector<unsigned int>::const_iterator count = kriging.GetLagCounts().begin();

		//assert(kriging.GetLagDistances().size() == kriging.GetLagSemivariances().size() && kriging.GetLagSemivariances().size() == kriging.GetLagCounts().size());

		//std::vector<double> linearVariograms = kriging.CalculateLinearModel();

		//std::vector<double>::iterator linearVariogram = linearVariograms.begin();

		//while (distance != kriging.GetLagDistances().end())
		//{
		//	std::cout << *distance++ << ", " << *semivariance++ << ", " << *count++ << ", " << *linearVariogram++ << std::endl;
		//}

		// kriging.SimpleKrige(Kriging::Spherical, 0.0, kriging.GetEstimatedSill(), 4000);

		kriging.OrdinaryKrige(_Kriging::Spherical, 0.0, kriging.GetEstimatedSill(), 4000, 16, 16, 0);

		/*
		std::cout << "Distance" << ", " << "Semivariance" << ", " << "Counts" << ", " <<  "Linear Variogram" << std::endl;

		std::vector<double>::const_iterator distance = kriging.GetLagDistances().begin();
		std::vector<double>::const_iterator semivariance = kriging.GetLagSemivariances().begin();
		std::vector<unsigned int>::const_iterator count = kriging.GetLagCounts().begin();

		assert(kriging.GetLagDistances().size() == kriging.GetLagSemivariances().size() && kriging.GetLagSemivariances().size() == kriging.GetLagCounts().size());

		std::vector<double> linearVariograms = kriging.CalculateLinearModel();

		std::vector<double>::iterator linearVariogram = linearVariograms.begin();

		while (distance != kriging.GetLagDistances().end())
		{
		std::cout << *distance++ << ", " << *semivariance++ << ", " << *count++ << ", " << *linearVariogram++ << std::endl;
		}
		*/

		// Note: Because of checked iterators, the debug build will run extremely slowly.  Consider turning
		//       this option off if debugging large datasets.

		// kriging.SimpleKrige(Kriging::Spherical, 0.0, kriging.GetEstimatedSill(), 4000);

		// kriging.OrdinaryKrige(Kriging::Spherical, 0.0, kriging.GetEstimatedSill(), 4000, 16, 16, 0.0);
	}

	return result;
}

void KrigeLargeSoilSampleSet()
{
	bool result = false;

	std::vector<DataPoint> dataPointP;
	std::vector<DataPoint> dataPointK;
	std::vector<DataPoint> dataPointCa;

	std::ifstream dataFile(".\\data\\LargeSoilSampleSet.csv", std::ifstream::in);

	if (dataFile.good())
	{
		result = true;

		std::string line;

		// Skip header row

		std::getline(dataFile, line);

		while (!dataFile.eof())
		{
			std::getline(dataFile, line);

			std::stringstream lineStream(line);

			std::string cell;

			unsigned int column = 0;

			double x = 0.0;
			double y = 0.0;

			double value = 0.0;

			while (std::getline(lineStream, cell, ','))
			{
				if (column == 0)
				{
					y = std::atof(cell.c_str());
				}
				else if (column == 1)
				{
					x = std::atof(cell.c_str());
				}
				else if (column == 3)
				{
					value = std::atof(cell.c_str());

					dataPointP.push_back(DataPoint(x, y, value));
				}
				else if (column == 4)
				{
					value = std::atof(cell.c_str());

					dataPointK.push_back(DataPoint(x, y, value));
				}
				else if (column == 5)
				{
					value = std::atof(cell.c_str());

					dataPointCa.push_back(DataPoint(x, y, value));
				}

				column++;
			}
		}
	}

	if (result)
	{
		_RasterContext rasterContext;

		_Kriging kriging(dataPointCa, rasterContext);

		kriging.Initialize();

		std::cout << "Distance" << ", " << "Semivariance" << ", " << "Counts" << ", " << "Linear Variogram" << std::endl;

		std::vector<double>::const_iterator distance = kriging.GetLagDistances().begin();
		std::vector<double>::const_iterator semivariance = kriging.GetLagSemivariances().begin();
		std::vector<unsigned int>::const_iterator count = kriging.GetLagCounts().begin();

		assert(kriging.GetLagDistances().size() == kriging.GetLagSemivariances().size() && kriging.GetLagSemivariances().size() == kriging.GetLagCounts().size());

		std::vector<double> linearVariograms = kriging.CalculateLinearModel();

		std::vector<double>::iterator linearVariogram = linearVariograms.begin();

		while (distance != kriging.GetLagDistances().end())
		{
			std::cout << *distance++ << ", " << *semivariance++ << ", " << *count++ << ", " << *linearVariogram++ << std::endl;
		}
	}

	return;
}


KrigingAlgo::RasterContext::RasterContext(double xOffsetSource, double yOffsetSource,
	double xScaleSource, double yScaleSource,
	unsigned int widthSource, unsigned int heightSource, bool yFlippedSource) {
	instance = new _RasterContext(xOffsetSource, yOffsetSource, xScaleSource, yScaleSource, widthSource, heightSource, yFlippedSource);

}

KrigingAlgo::RasterContext::~RasterContext() {
	if (instance) {
		_RasterContext* ptr = static_cast<_RasterContext*> (instance);
		delete ptr;
	}
}

KrigingAlgo::Kriging::Kriging(array<double>^ x, array<double>^ y, array<double>^ z, RasterContext^ rasterContext)
{
	if (x->Length == y->Length &&
		x->Length == z->Length &&
		y->Length == z->Length) {
		std::vector<DataPoint> dataPoint;
		for (int i = 0; i < x->Length; i++)
		{
			dataPoint.push_back(DataPoint(x[i], y[i], z[i]));
		}
		_RasterContext* ptr = static_cast<_RasterContext*>(rasterContext->instance);
		instance = new _Kriging(dataPoint, *ptr);
	}
	else {
		throw gcnew System::Exception("xyz数组长度不一致");
	}
}

KrigingAlgo::Kriging::~Kriging()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		delete ptr;
	}
}

void KrigingAlgo::Kriging::Initialize()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		ptr->Initialize();
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

void KrigingAlgo::Kriging::CalculateExperimentalVariogram(double lag, double lagTolerance)
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		ptr->CalculateExperimentalVariogram(lag, lagTolerance);
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

array<double>^ KrigingAlgo::Kriging::CalculateLinearModel()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		auto items = ptr->CalculateLinearModel();
		array<double>^ arr = gcnew array<double>(items.size());
		for (int i = 0; i < items.size(); i++)
		{
			arr[i] = items[i];
		}
		return arr;
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

array<double>^ KrigingAlgo::Kriging::CalculateLinearModelWithoutIntercept(double nugget)
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		auto items = ptr->CalculateLinearModelWithoutIntercept(nugget);
		array<double>^ arr = gcnew array<double>(items.size());
		for (int i = 0; i < items.size(); i++)
		{
			arr[i] = items[i];
		}
		return arr;
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

array<double>^ KrigingAlgo::Kriging::CalculateSphericalModel(double nugget, double sill, double range)
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		auto items = ptr->CalculateSphericalModel(nugget, sill, range);
		array<double>^ arr = gcnew array<double>(items.size());
		for (int i = 0; i < items.size(); i++)
		{
			arr[i] = items[i];
		}
		return arr;
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

array<double>^ KrigingAlgo::Kriging::CalculateExponentialModel(double nugget, double sill, double range)
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		auto items = ptr->CalculateExponentialModel(nugget, sill, range);
		array<double>^ arr = gcnew array<double>(items.size());
		for (int i = 0; i < items.size(); i++)
		{
			arr[i] = items[i];
		}
		return arr;
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

array<double>^ KrigingAlgo::Kriging::CalculateGaussianModel(double nugget, double sill, double range)
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		auto items = ptr->CalculateGaussianModel(nugget, sill, range);
		array<double>^ arr = gcnew array<double>(items.size());
		for (int i = 0; i < items.size(); i++)
		{
			arr[i] = items[i];
		}
		return arr;
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

double KrigingAlgo::Kriging::GetEstimatedNugget()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		return ptr->GetEstimatedNugget();
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

double KrigingAlgo::Kriging::GetEstimatedRange()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		return ptr->GetEstimatedRange();
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

double KrigingAlgo::Kriging::GetEstimatedSill()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		return ptr->GetEstimatedSill();
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

double KrigingAlgo::Kriging::GetDefaultLag()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		return ptr->GetDefaultLag();
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

double KrigingAlgo::Kriging::GetDefaultLagTolerance()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		return ptr->GetDefaultLagTolerance();
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

array<double>^ KrigingAlgo::Kriging::GetLagDistances()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		auto items = ptr->GetLagDistances();
		array<double>^ arr = gcnew array<double>(items.size());
		for (int i = 0; i < items.size(); i++)
		{
			arr[i] = items[i];
		}
		return arr;
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

array<double>^ KrigingAlgo::Kriging::GetLagSemivariances()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		auto items = ptr->GetLagSemivariances();
		array<double>^ arr = gcnew array<double>(items.size());
		for (int i = 0; i < items.size(); i++)
		{
			arr[i] = items[i];
		}
		return arr;
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

array<unsigned int>^ KrigingAlgo::Kriging::GetLagCounts()
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		auto items = ptr->GetLagCounts();
		array<unsigned int>^ arr = gcnew array<unsigned int>(items.size());
		for (int i = 0; i < items.size(); i++)
		{
			arr[i] = items[i];
		}
		return arr;
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

void KrigingAlgo::Kriging::SimpleKrige(Model model, array<float>^ buffer, double nugget, double sill, double range)
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		_Kriging::_Model m = _Kriging::_Model::Linear;
		switch (model)
		{
		case KrigingAlgo::Kriging::Model::Linear:
			m = _Kriging::_Model::Linear;
			break;
		case KrigingAlgo::Kriging::Model::LinearWithoutIntercept:
			m = _Kriging::_Model::LinearWithoutIntercept;
			break;
		case KrigingAlgo::Kriging::Model::Spherical:
			m = _Kriging::_Model::Spherical;
			break;
		case KrigingAlgo::Kriging::Model::Exponential:
			m = _Kriging::_Model::Exponential;
			break;
		case KrigingAlgo::Kriging::Model::Gaussian:
			m = _Kriging::_Model::Gaussian;
			break;
		default:
			break;
		}
		GCHandle handle = GCHandle::Alloc(buffer, GCHandleType::Pinned);
		IntPtr bufferPtr = handle.AddrOfPinnedObject();
		
		ptr->SimpleKrige(m, bufferPtr.ToPointer(), nugget, sill, range);
		handle.Free();
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}

void KrigingAlgo::Kriging::OrdinaryKrige(Model model, array<float>^ buffer, double nugget, double sill, double range, unsigned int minPoints, unsigned int maxPoints, double maxDistance)
{
	if (instance) {
		_Kriging* ptr = static_cast<_Kriging*>(instance);
		_Kriging::_Model m = _Kriging::_Model::Linear;
		switch (model)
		{
		case KrigingAlgo::Kriging::Model::Linear:
			m = _Kriging::_Model::Linear;
			break;
		case KrigingAlgo::Kriging::Model::LinearWithoutIntercept:
			m = _Kriging::_Model::LinearWithoutIntercept;
			break;
		case KrigingAlgo::Kriging::Model::Spherical:
			m = _Kriging::_Model::Spherical;
			break;
		case KrigingAlgo::Kriging::Model::Exponential:
			m = _Kriging::_Model::Exponential;
			break;
		case KrigingAlgo::Kriging::Model::Gaussian:
			m = _Kriging::_Model::Gaussian;
			break;
		default:
			break;
		}
		GCHandle handle = GCHandle::Alloc(buffer, GCHandleType::Pinned);
		IntPtr bufferPtr = handle.AddrOfPinnedObject();

		ptr->OrdinaryKrige(m, bufferPtr.ToPointer(), nugget, sill, range, minPoints, maxPoints, maxDistance);
		handle.Free();
	}
	else
	{
		throw gcnew System::NullReferenceException("克里金对象指针为空");
	}
}


