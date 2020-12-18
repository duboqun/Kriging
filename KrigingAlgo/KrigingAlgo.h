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
		double xOffsetSource; double yOffsetSource;

		double xScaleSource; double yScaleSource;
		int widthSource; int heightSource;
	};
	public enum class Model
	{
		Linear,
		LinearWithoutIntercept,
		Spherical,
		Exponential,
		Gaussian
	};

	public ref class Kriging
	{
	private:
		void* instance = nullptr;
	public:


		Kriging(cli::array<double>^ x, cli::array<double>^ y, cli::array<double>^ z, RasterContext^ rasterContext);
		~Kriging();
		void Initialize();

		void CalculateExperimentalVariogram(double lag, double lagTolerance);
		cli::array<double>^ CalculateLinearModel();
		cli::array<double>^ CalculateLinearModelWithoutIntercept(double nugget);
		cli::array<double>^ CalculateSphericalModel(double nugget, double sill, double range);
		cli::array<double>^ CalculateExponentialModel(double nugget, double sill, double range);
		cli::array<double>^ CalculateGaussianModel(double nugget, double sill, double range);

		
		double GetEstimatedNugget();
		double GetEstimatedRange();
		double GetEstimatedSill();
		double GetDefaultLag();
		double GetDefaultLagTolerance();
		cli::array<double>^ GetLagDistances();
		cli::array<double>^ GetLagSemivariances();
		cli::array<unsigned int>^ GetLagCounts();


		void SimpleKrige(Model model,cli::array<float>^ buffer, double nugget, double sill, double range);
		void OrdinaryKrige(Model model, cli::array<float>^ buffer, double nugget, double sill, double range,
			unsigned int minPoints, unsigned int maxPoints, double maxDistance);

	};
}
