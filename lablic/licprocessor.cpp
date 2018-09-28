/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <lablic/licprocessor.h>
#include <lablic/integrator.h>
#include <lablic/interpolator.h>
#include <inviwo/core/datastructures/volume/volumeram.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
	: Processor()
	, volumeIn_("volIn")
	, noiseTexIn_("noiseTexIn")
	, licOut_("licOut")
	, kernelLength("klength", "Length of kernel", 0, 0, 10)

// TODO: Register additional properties

{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties

    // TODO: Register additional properties
	addProperty(kernelLength);

}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    auto vol = volumeIn_.getData();
    vectorFieldDims_ = vol->getDimensions();
    auto vr = vol->getRepresentation<VolumeRAM>();

    // An accessible form of on image is retrieved analogous to a volume
    auto tex = noiseTexIn_.getData();
    texDims_ = tex->getDimensions();
    auto tr = tex->getRepresentation<ImageRAM>();

    // Prepare the output, it has the same dimensions and datatype as the output
    // and an editable version is retrieved analogous to a volume
    auto outImage = tex->clone();
    auto outLayer = outImage->getColorLayer();
    auto lr = outLayer->getEditableRepresentation<LayerRAM>();

    // To access the image at a floating point position, you can call
    //      Interpolator::sampleFromGrayscaleImage(tr, somePos)

    // TODO: Implement LIC and FastLIC
	//std::vector<std::vector<int>> visited(texDims_.x, std::vector<int>(texDims_.y, 0));
	for (auto j = 0; j < texDims_.y-1; j++) {
		for (auto i = 0; i < texDims_.x-1; i++) {
			int val = applyKernel(vec2(i, j), vol.get(), tr);
			lr->setFromDVec4(size2_t(i, j), dvec4(val, val, val, 255));
		}
	}


    // This code instead sets all pixels to the same gray value
	/*
    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, 255));

    for (auto j = 0; j < texDims_.y; j++) {
        for (auto i = 0; i < texDims_.x; i++) {
            int val = int(licTexture[i][j]);
            lr->setFromDVec4(size2_t(i, j), dvec4(val, val, val, 255));
        }
    }
	*/

    licOut_.setData(outImage);
}

std::vector<double> LICProcessor::kernel() {
	float weigth = 1. / (kernelLength.get() * 2 + 1);
	return std::vector<double>(kernelLength.get(), weigth);
}

int LICProcessor::applyKernel(const vec2& pos, const Volume* vr, const ImageRAM* tr) {
	//go kernel length forward and backward
	double sum = Interpolator::sampleFromGrayscaleImage(tr, pos);
	vec2 newPos = textureToField(pos);
	vec2 dir;
	int div = 1;
	for (int i = 0; i < kernelLength.get(); i++) {
		dir = Integrator::RK4(vr, newPos, 1., 1.);
		if (!(Integrator::vecLength(dir) < 0.001)) {
			div += 1;
			newPos += dir;
			vec2 texPos = fieldToTexture(newPos);
			sum += Interpolator::sampleFromGrayscaleImage(tr, texPos);
		}
		else
			break;
	}
	newPos = textureToField(pos);
	for (int i = 0; i < kernelLength.get(); i++) {
		dir = Integrator::RK4(vr, newPos, 1., -1.);
		if (!(Integrator::vecLength(dir) < 0.001)) {
			div += 1;
			newPos += dir;
			vec2 texPos = fieldToTexture(newPos);
			sum += Interpolator::sampleFromGrayscaleImage(tr, texPos);
		}
		else 
			break;	
	}
	return (int)(sum / div);
}

vec2 LICProcessor::fieldToTexture(const vec2& pos) {
	double xRatio = (double)vectorFieldDims_.x / (double)texDims_.x;
	double yRatio = (double)vectorFieldDims_.y / (double)texDims_.y;
	return vec2((int)clip(pos.x / xRatio,0,texDims_.x-2), (int)clip(pos.y / yRatio,0,texDims_.y-2));
}

vec2 LICProcessor::textureToField(const vec2& pos) {
	double xRatio = (double)texDims_.x / (double)vectorFieldDims_.x;
	double yRatio = (double)texDims_.y / (double)vectorFieldDims_.y;
	return vec2(pos.x / xRatio, pos.y / yRatio);
}

double LICProcessor::clip(double n, double lower, double upper) {
	return std::max(lower, std::min(n, upper));
}


}  // namespace inviwo
