#ifndef CRPROPA_AMRMAGNETICFIELD_H
#define CRPROPA_AMRMAGNETICFIELD_H

#ifdef CRPROPA_HAVE_SAGA

#include <iostream> 
#include <string>
#include <cstdio>

#ifdef _OPENMP
    #include "omp.h"
#endif

#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Vector3.h"

#include "saga/LocalProperties.h"
#include "saga/AMRgrid.h"
#include "saga/MagneticField.h"
#include "saga/Referenced.h"



namespace crpropa {

/**
 @class AMRMagneticField
 @brief Wrapper for saga::MagneticField
 */
class AMRMagneticField: public MagneticField {

private:
	saga::ref_ptr<saga::MagneticField> field;
    double cfLength;
    double cfDensity;
    double cfMagneticField;

public:        
  
    AMRMagneticField(saga::ref_ptr<saga::MagneticField> field_, double convLength, double convDensity, double convMagneticField)            
    {
        field = field_;
        cfLength = convLength;
        cfDensity = convDensity;
        cfMagneticField = convMagneticField;
    }

    Vector3d getField(const Vector3d &position) const {

        double x = position.x/cfLength;
        double y = position.y/cfLength;
        double z = position.z/cfLength;

        std::vector<double> b;
        #ifdef _OPENMP
            #pragma omp critical 
            {    
		        b = field->getField(x, y, z);
            }
        #else 
            b = field->getField(x, y, z);
        #endif

        for(int i=0; i<3; i++)
            b[i]*=cfMagneticField;
        //std::cout << x*cfLength/Mpc << "  " << y*cfLength/Mpc << "  " << z*cfLength/Mpc << "  ||   " << b[0] << " " << b[1] << "  " << b[0] << std::endl;

		return Vector3d(b[0], b[1], b[2]) * tesla;   
    }

};

} // namespace crpropa

#endif // CRPROPA_HAVE_SAGA
#endif // CRPROPA_AMRMAGNETICFIELD_H
