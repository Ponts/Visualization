/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <ctime>
#include <limits>
#include <lablic/lablicmoduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/ports/imageport.h>
#include <inviwo/core/datastructures/image/imageram.h>
#include <inviwo/core/properties/boolproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>



namespace inviwo {

/** \docpage{org.inviwo.LICProcessor, LICProcessor}
    ![](org.inviwo.LICProcessor.png?classIdentifier=org.inviwo.LICProcessor)

    Line Integral Convolution with a box kernel.
    
    ### Inports
      * __vectorField__ 2-dimensional vector field (with vectors of
      two components thus two values within each voxel) represented
      by a 3-dimensional volume.
      This processor deals with 2-dimensional data only, therefore it is assumed
      the z-dimension will have size 1 otherwise the 0th slice of the volume
      will be processed.
      * __texture__ Texture to be convolved along the streamlines.
    
    ### Outports
      * __image__ The image resulting from smearing the given texture
      the streamlines of the given vector field.
*/
class IVW_MODULE_LABLIC_API LICProcessor : public Processor {
    // Friends
// Types
public:
// Construction / Deconstruction
public:
    LICProcessor();
    virtual ~LICProcessor() = default;

// Methods
public:
    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;

protected:
    // Our main computation function
    virtual void process() override;
	std::vector<double> kernel();
	int applyKernel(const vec2& pos, const Volume* vr, const ImageRAM* tr);
	void fastKernel(const vec2& pos, const Volume* vr, const ImageRAM* tr, std::vector<std::vector<int>>& visited, LayerRAM* lr);
	vec2 fieldToTexture(const vec2& pos);
	vec2 textureToField(const vec2& pos);
	double clip(double n, double lower, double upper);
	void enhance(LayerRAM* lr);
	void colorise(LayerRAM* lr, const Volume* vr);
	vec4 getColor(const double mag);

    // (TODO: Helper functions can be defined here and then implemented in the .cpp)
    // e.g. something like a function for standardLIC, fastLIC, autoContrast, ...

// Ports
public:
    // Input vector field
    VolumeInport volumeIn_;

    // Input texture
    ImageInport noiseTexIn_;

    // Output image
    ImageOutport licOut_;

// Properties
public:

// TODO: Declare properties
// IntProperty prop1;
// BoolProperty prop2;
	IntProperty kernelLength;
	BoolProperty improveContrast;
	BoolProperty fastLic;
	BoolProperty debug;
	BoolProperty colorIt;

// Attributes
private:
    size3_t vectorFieldDims_;
    size2_t texDims_;
	double fieldToTexXrat;
	double fieldToTexYrat;

	double texToFieldXrat;
	double texToFieldYrat;

	double stepSize;

	


};

}  // namespace inviwo
