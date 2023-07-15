#ifndef LPSI4_H //Usual macro guard to prevent multiple inclusion
#define LPSI4_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2020

#include "CosmoInterface/cosmointerface.h"

//Include cosmointerface to have access to all of the library.

namespace TempLat
{
    /////////
    // Model name and number of fields
    /////////

    // In the following class, we define the defining parameters of your model:
    // number of fields of each species and the type of tinteractions.

    struct ModelPars : public TempLat::DefaultModelPars {
    	static constexpr size_t NScalars = 2;
        static constexpr size_t NPotTerms = 3;
    };

  #define MODELNAME lpsi4
  // Here we define the name of the model. This should match the name of your file.

  template<class R>
  using Model = MakeModel(R, ModelPars);
  // In this line, we define an appropriate generic model, with the correct
  // number of fields, ready to be customized.
  // If you are curious about what this is doing, the macro is defined in
  // the "CosmoInterface/abstractmodel.h" file.

  class MODELNAME : public Model<MODELNAME>
  // Declaration of our model. It inherits from the generic model defined above.
  {
 //...
private:

  double big_m, m, lambda, gamma, q, beta_sq;

  public:

    MODELNAME(ParameterParser& parser, RunParameters<double>& runPar, std::shared_ptr<MemoryToolBox> toolBox): //Constructor of our model.
    Model<MODELNAME>(parser,runPar.getLatParams(), toolBox, runPar.dt, STRINGIFY(MODELLABEL)) //MODELLABEL is defined in the cmake.
    {

      /////////
      // Independent parameters of the model (read from parameters file)
      /////////
      big_m = parser.get<double>("big_m");
      m = parser.get<double>("m");
      lambda = parser.get<double>("lambda");
      gamma = parser.get<double>("gamma");

      beta_sq = pow<2>(m/big_m)*1/(2*sqrt(lambda));
      q = gamma/lambda;

      std::cout << "beta_sq = " << beta_sq << " q = " << q << std::endl;
        /////////
        // Initial homogeneous components of the fields
        // (read from parameters file, or specified here if not)
        /////////

        fldS0 = parser.get<double, 2>("initial_amplitudes");
        piS0 = parser.get<double, 2>("initial_momenta", {0, 0});

        /////////
        // Rescaling for program variables
        /////////

        alpha = 0;
        fStar = sqrt(2/sqrt(lambda))*big_m;
        omegaStar = sqrt(2*sqrt(lambda))*big_m;

        setInitialPotentialAndMassesFromPotential();
    }

   /////////
   // Program potential (add as many functions as terms are in the potential)
   /////////

    auto potentialTerms(Tag<0>)
    {
        return 0.5 *beta_sq * pow<2>(fldS(0_c));
    }
    auto potentialTerms(Tag<1>) // Interaction energy
    {
        return 0.25  * pow<4>(fldS(1_c)) - 0.5 * pow<2>(fldS(1_c));
    }
    auto potentialTerms(Tag<2>)
    {
        return 0.5 * q * pow<2>(fldS(0_c) * fldS(1_c)) + 0.25 ;
    }
   /////////
   // Derivatives of the program potential with respect fields
   // (add one function for each field).
   /////////

    auto potDeriv(Tag<0>) 
    {
      return  beta_sq*fldS(0_c) + q * pow<2>(fldS(1_c)) * fldS(0_c) ;
    }

    auto potDeriv(Tag<1>)
    {
      return pow<3>(fldS(1_c)) - fldS(1_c) + q * pow<2>(fldS(0_c)) * fldS(1_c) ;
    }
	
	
	
    /////////
   //  Second derivatives of the program potential with respect fields
   // (add one function for each field)
   /////////

    auto potDeriv2(Tag<0>)
    {
      return  beta_sq + q* pow<2>(fldS(1_c)) ;
    }

    auto potDeriv2(Tag<1>)
    {
      return  3 * pow<2>(fldS(1_c)) - 1 + q * pow<2>(fldS(0_c)) ;
    }
		
    };
}

#endif //LPSI4_H
