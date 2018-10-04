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
	, kernelLength("klength", "Length of kernel", 0, 0, 100)
	, improveContrast("icontr", "Improve Contrast", false)
	, fastLic("fast", "Fast LIC", true)
	, debug("go", "go", false)
	, colorIt("color", "Colorize", false)

// TODO: Register additional properties

{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties

    // TODO: Register additional properties
	addProperty(kernelLength);
	addProperty(improveContrast);
	addProperty(fastLic);
	addProperty(debug);
	addProperty(colorIt);

}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

	if (!debug) {
		return;
	}

    auto vol = volumeIn_.getData();
    vectorFieldDims_ = vol->getDimensions();
    //auto vr = vol->getRepresentation<VolumeRAM>();

    // An accessible form of on image is retrieved analogous to a volume
    auto tex = noiseTexIn_.getData();
    texDims_ = tex->getDimensions();
    auto tr = tex->getRepresentation<ImageRAM>();

    // Prepare the output, it has the same dimensions and datatype as the output
    // and an editable version is retrieved analogous to a volume
    auto outImage = tex->clone();
    auto outLayer = outImage->getColorLayer();
    auto lr = outLayer->getEditableRepresentation<LayerRAM>();
	fieldToTexXrat = (double)(vectorFieldDims_.x - 1) / (double)(texDims_.x - 1);
	fieldToTexYrat = (double)(vectorFieldDims_.y - 1) / (double)(texDims_.y - 1);

	texToFieldXrat = (double)(texDims_.x - 1) / (double)(vectorFieldDims_.x - 1);
	texToFieldYrat = (double)(texDims_.y - 1) / (double)(vectorFieldDims_.y - 1);

	stepSize = std::min(fieldToTexXrat, fieldToTexYrat);

    // To access the image at a floating point position, you can call
    //      Interpolator::sampleFromGrayscaleImage(tr, somePos)

    // TODO: Implement LIC and FastLIC
	int entered = 0;
	auto volr = vol.get();
	std::time_t start = std::time(nullptr);
	if (fastLic) {
		std::vector<std::vector<int>> visited(texDims_.x, std::vector<int>(texDims_.y, 0));
		int xCellSize = texDims_.x / 6;
		int yCellSize = texDims_.y / 6;
		for (auto j = 0; j < yCellSize; j++) {
			for (auto i = 0; i < xCellSize; i++) {
				for (auto blockY = 0; blockY < texDims_.y; blockY += yCellSize) {
					for (auto blockX = 0; blockX < texDims_.x; blockX += xCellSize) {
						int x = blockX + i;
						int y = blockY + j;
						if (texDims_.x <= x || texDims_.y <= y) continue;
						if (visited[x][y] == 0) {
							fastKernel(vec2(x, y), volr, tr, visited, lr);
							entered++;
						}
					}
				}
				
			}
		}
	}
	else {
		for (auto j = 0; j < texDims_.y; j++) {
			for (auto i = 0; i < texDims_.x; i++) {
				int val = applyKernel(vec2(i, j), volr, tr);
				lr->setFromDVec4(size2_t(i, j), dvec4(val, val, val, 255));
				entered++;
			}
		}
	}
	std::time_t stop = std::time(nullptr);
	LogProcessorInfo(entered);
	LogProcessorInfo(stop - start << "s run time");
	


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
	if (improveContrast)
		enhance(lr);
	
	if (colorIt)
		colorise(lr, vol.get());
    licOut_.setData(outImage);
}

void LICProcessor::colorise(LayerRAM* lr, const Volume* vr) {
	//find max magnitude
	double maxMag = -1.;
	for (auto j = 0; j < vectorFieldDims_.y; j++) {
		for (auto i = 0; i < vectorFieldDims_.x; i++) {
			double mag = Integrator::vecLength(Interpolator::sampleFromField(vr,vec2(i,j)));
			if (maxMag < mag)
				maxMag = mag;
		}
	}
	for (auto j = 0; j < texDims_.y; j++) {
		for (auto i = 0; i < texDims_.x; i++) {
			double licFactor = lr->getAsDouble(size2_t(i,j)) / 255;
			double normMag = Integrator::vecLength(Interpolator::sampleFromField(vr, textureToField(vec2(i, j)))) / maxMag;
			dvec4 color = getColor(normMag);
			for (int c = 0; c < 3; c++)
				color[c] = color[c] * licFactor;
			lr->setFromDVec4(size2_t(i, j), color);
		}
	}
}

vec4 LICProcessor::getColor(const double mag) {
	if (mag <= 1e-20)
		return vec4(0, 0, 0, 255);
	return vec4((int)(mag*255), 0, (int)((1-mag)*255), 255);
}

void LICProcessor::enhance(LayerRAM* lr) {
	double mu = 0.;
	double P = 0.;
	double n = 0.;
	for (auto j = 0; j < texDims_.y; j++) {
		for (auto i = 0; i < texDims_.x; i++) {
			double p = lr->getAsDouble(size2_t(i, j));
			if (p > 0.001) {
				mu += p;
				n+=1.;
				P += pow(p, 2);
			}
		}
	}
	mu /= n;
	double sigma = sqrt((P - n*pow(mu, 2)) / (n - 1));
	double desMean = 255 * 0.5;
	double desVar = 0.2 * 255;
	double f = desVar / sigma;
	for (auto j = 0; j < texDims_.y; j++) {
		for (auto i = 0; i < texDims_.x; i++) {
			double p = lr->getAsDouble(size2_t(i, j));
			double pi = desMean + (f*(p - mu));
			lr->setFromDVec4(size2_t(i, j), dvec4(pi, pi, pi, 255));
		}
	}
	return;
}

std::vector<double> LICProcessor::kernel() {
	float weigth = 1. / (kernelLength.get() * 2 + 1);
	return std::vector<double>(kernelLength.get(), weigth);
}

int LICProcessor::applyKernel(const vec2& pos, const Volume* vr, const ImageRAM* tr) {
	//go kernel length forward and backward
	double sum = tr->readPixel(size2_t(pos[0], pos[1]), LayerType::Color)[0];
	vec2 newPos = textureToField(pos);
	vec2 dir;
	int div = 1;
	for (int i = 0; i < kernelLength.get(); i++) {
		dir = Integrator::RK4(vr, newPos, stepSize, 1.);
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
		dir = Integrator::RK4(vr, newPos, stepSize, -1.);
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

void LICProcessor::fastKernel(const vec2& pos, const Volume* vr, const ImageRAM* tr, std::vector<std::vector<int>>& visited, LayerRAM* lr) {
	double sum = tr->readPixel(size2_t(pos[0], pos[1]), LayerType::Color)[0];
	std::queue<vec2> forward;
	std::queue<vec2> backward;
	std::vector<double> values(kernelLength*2+1,0);
	int div = 1;
	vec2 newPos = textureToField(vec2(pos.x,pos.y)), forwardPos, backPos, dir;
	values[kernelLength] = sum;
	//forward
	for (int i = 0; i < kernelLength.get(); i++) {
		dir = Integrator::RK4(vr, newPos, stepSize, 1.);
		if (!(Integrator::vecLength(dir) < 0.001)) {
			div += 1;
			newPos += dir;
			forward.push(newPos);
			vec2 texPos = fieldToTexture(newPos);
			double value = Interpolator::sampleFromGrayscaleImage(tr, texPos);
			values[kernelLength + i + 1] = value;
			sum += value;
		}
		else 
			break;
	}
	newPos = textureToField(pos);
	for (int i = 0; i < kernelLength.get(); i++) {
		dir = Integrator::RK4(vr, newPos, stepSize, -1.);
		if (!(Integrator::vecLength(dir) < 0.001)) {
			div += 1;
			newPos += dir;
			backward.push(newPos);
			vec2 texPos = fieldToTexture(newPos);
			double value = Interpolator::sampleFromGrayscaleImage(tr, texPos);
			values[kernelLength - i - 1] = value;
			sum += value;
		}
		else
			break;
	}
	double startSum = sum;
	double startDiv = div;
	std::vector<double> startValues = values;
	std::queue<double> forwardValues;
	for (const auto& e : values) forwardValues.push(e);
	visited[pos.x][pos.y] += 1; //visited
	double color = sum / (double)div;
	lr->setFromDVec4(size2_t(pos.x, pos.y), vec4(color, color, color, 255));
	if (forward.size() != 0) {
		dir = Integrator::RK4(vr, forward.back(), stepSize, 1.);
		while (!(Integrator::vecLength(dir) < 0.001) ) {
			vec2 colorPos = fieldToTexture(forward.front());
			if (visited[(int)colorPos.x][(int)colorPos.y] >= 10) break;
			newPos = forward.back() + dir;
			forward.push(newPos);
			double value = Interpolator::sampleFromGrayscaleImage(tr, fieldToTexture(newPos));
			sum += value;
			sum -= forwardValues.front();
			forwardValues.pop();
			forwardValues.push(value);
			div = std::min(div + 1, kernelLength * 2 + 1);
			if (visited[(int)colorPos.x][(int)colorPos.y] < 1) {
				double color = sum / (double)div;
				lr->setFromDVec4(size2_t((int)colorPos.x, (int)colorPos.y), vec4(color, color, color, 255));
			}
			visited[(int)colorPos.x][(int)colorPos.y] += 1;
			forward.pop();
			dir = Integrator::RK4(vr, forward.back(), stepSize, 1.);
		}
		//middle of kernel hit visited or edge
		for (const auto& e = forward.front(); !forward.empty(); forward.pop()) {
			vec2 imagePos = fieldToTexture(e);
			if (visited[(int)imagePos.x][(int)imagePos.y] >= 10) break;
			if (visited[(int)imagePos.x][(int)imagePos.y] < 1) {
				double color = sum / (double)div;
				lr->setFromDVec4(size2_t((int)imagePos.x, (int)imagePos.y), vec4(color, color, color, 255));
			}
			visited[(int)imagePos.x][(int)imagePos.y] += 1;

			sum -= forwardValues.front();
			forwardValues.pop();
			div -= 1;
		}
	}
	sum = startSum;
	div = startDiv;
	values = startValues;
	std::queue<double> backwardValues;
	for (int i = values.size() - 1; i >= 0; i-- ) backwardValues.push(values[i]);
	if (backward.size() != 0) {
		dir = Integrator::RK4(vr, backward.back(), stepSize, -1.);
		while (!(Integrator::vecLength(dir) < 0.001)) {
			vec2 colorPos = fieldToTexture(backward.front());
			if (visited[(int)colorPos.x][(int)colorPos.y] >= 10) break;
			newPos = backward.back() + dir;
			backward.push(newPos);
			double value = Interpolator::sampleFromGrayscaleImage(tr, fieldToTexture(newPos));
			sum += value;
			sum -= backwardValues.front();
			backwardValues.pop();
			backwardValues.push(value);
			div = std::min(div + 1, kernelLength * 2 + 1);
			if (visited[(int)colorPos.x][(int)colorPos.y] < 1) {
				double color = sum / (double)div;
				lr->setFromDVec4(size2_t((int)colorPos.x, (int)colorPos.y), vec4(color, color, color, 255));
			}
			visited[(int)colorPos.x][(int)colorPos.y] += 1;
			backward.pop();
			dir = Integrator::RK4(vr, backward.back(), stepSize, -1.);
		}
		//middle of kernel hit visited or edge
		for (const auto& e = backward.front(); !backward.empty(); backward.pop()) {
			vec2 imagePos = fieldToTexture(e);
			if (visited[(int)imagePos.x][(int)imagePos.y] >= 10) break;
			if (visited[(int)imagePos.x][(int)imagePos.y] < 1) {
				double color = sum / (double)div;
				lr->setFromDVec4(size2_t((int)imagePos.x, (int)imagePos.y), vec4(color, color, color, 255));
			}
			visited[(int)imagePos.x][(int)imagePos.y] += 1;
			sum -= backwardValues.front();
			backwardValues.pop();
			div -= 1;
		}
	}

}

vec2 LICProcessor::fieldToTexture(const vec2& pos) {
	return vec2(clip(pos.x / fieldToTexXrat,0,texDims_.x-1), clip(pos.y / fieldToTexYrat,0,texDims_.y-1));
}

vec2 LICProcessor::textureToField(const vec2& pos) {
	return vec2(clip(pos.x / texToFieldXrat,0, vectorFieldDims_.x-1), clip(pos.y / texToFieldYrat,0, vectorFieldDims_.y-1));
}

double LICProcessor::clip(double n, double lower, double upper) {
	return std::max(lower, std::min(n, upper - 0.01));
}

//std::numeric_limits<double>::min()
}  // namespace inviwo
