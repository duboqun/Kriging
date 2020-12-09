// KrigingAlgo.h

#pragma once

using namespace System;

namespace KrigingAlgo {

	public ref class RasterContext {
	internal:
		void* instance = nullptr;
	public:
		RasterContext(double xOffsetSource, double yOffsetSource,
			double xScaleSource, double yScaleSource,
			unsigned int widthSource, unsigned int heightSource, bool yFlippedSource);

		~RasterContext();
	};

	public ref class Kriging
	{
	private:
		void* instance = nullptr;
	public:
		enum class Model
		{
			Linear,
			LinearWithoutIntercept,
			Spherical,
			Exponential,
			Gaussian
		};

		Kriging(array<double>^ x, array<double>^ y, array<double>^ z, RasterContext^ rasterContext);
		~Kriging();
		void Initialize();

		void CalculateExperimentalVariogram(double lag, double lagTolerance);
		array<double>^ CalculateLinearModel();
		array<double>^ CalculateLinearModelWithoutIntercept(double nugget);
		array<double>^ CalculateSphericalModel(double nugget, double sill, double range);
		array<double>^ CalculateExponentialModel(double nugget, double sill, double range);
		array<double>^ CalculateGaussianModel(double nugget, double sill, double range);

		
		double GetEstimatedNugget();
		double GetEstimatedRange();
		double GetEstimatedSill();
		double GetDefaultLag();
		double GetDefaultLagTolerance();
		array<double>^ GetLagDistances();
		array<double>^ GetLagSemivariances();
		array<unsigned int>^ GetLagCounts();


		void SimpleKrige(Model model,array<float>^ buffer, double nugget, double sill, double range);
		void OrdinaryKrige(Model model, array<float>^ buffer, double nugget, double sill, double range,
			unsigned int minPoints, unsigned int maxPoints, double maxDistance);

	};
}
