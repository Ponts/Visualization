/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <inviwo/core/datastructures/image/imageram.h>
#include <lablic/lablicmoduledefine.h>

namespace inviwo
{

class IVW_MODULE_LABLIC_API Interpolator {
    //Construction / Deconstruction
  public:
      Interpolator() {}
      virtual ~Interpolator() = default;

    //Methods
  public:
      static vec2 sampleFromField(const Volume* vol, const vec2& position);
      static double sampleFromGrayscaleImage(const ImageRAM* img, const vec2& position);
};

}// namespace inviwo
