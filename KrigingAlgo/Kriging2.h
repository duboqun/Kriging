#pragma once
#include "KrigingAlgo.h"
//Matrix Handling
#include "Matrix.h"
#include "Vector3.h"
//Local Files
#include "Utils.h"
#include "Variogram.h"
using namespace System;

namespace KrigingAlgo {
	public ref class Kriging2
	{
	private:
		vector3* vecArray = nullptr;
		int vecLength = 0;
		RasterContext^ rasterContext;
	public:

		Kriging2(cli::array<double>^ xArray, cli::array<double>^ yArray, cli::array<double>^ zArray, RasterContext^ rasterContext);
		~Kriging2();
		//void Initialize();


		void SimpleKrige(Model model, cli::array<float>^ buffer, double nugget, double sill, double range);
	};

}