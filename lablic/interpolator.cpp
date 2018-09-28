/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <lablic/interpolator.h>
#include <inviwo/core/datastructures/volume/volume.h>

namespace inviwo {

vec2 Interpolator::sampleFromField(const Volume* vol, const vec2& position) {

    auto vr = vol->getRepresentation<VolumeRAM>();
    auto dims = vr->getDimensions();

    auto base = vol->getBasis();
    vec2 cellSize(base[0][0] / dims[0], base[1][1] / dims[1]);

    // Sampled outside the domain!
    if (position[0] < 0 || position[0] > dims[0] - 1 || position[1] < 0 ||
        position[1] > dims[1] - 1) {
        return vec2(0, 0);
    }
    size2_t p((size_t)position[0], (size_t)position[1]);

    // Leads to accessing only inside the volume
    // Coefficients computation takes care of using the correct values
    for (int d = 0; d < 2; ++d) p[d] = std::min(p[d], dims[d] - 2);

    const auto f00 = vr->getAsDVec2(size3_t(p[0], p[1], 0));
    const auto f10 = vr->getAsDVec2(size3_t(p[0] + 1, p[1], 0));
    const auto f01 = vr->getAsDVec2(size3_t(p[0], p[1] + 1, 0));
    const auto f11 = vr->getAsDVec2(size3_t(p[0] + 1, p[1] + 1, 0));

    const float x = position[0] - p[0];
    const float y = position[1] - p[1];

    vec2 f;

    for (int i = 0; i < 2; i++) {
        f[i] = f00[i] * (1 - x) * (1 - y) + f01[i] * (1 - x) * y + f10[i] * x * (1 - y) +
               f11[i] * x * y;

        // Bring vector back to grid space.
        f[i] /= cellSize[i];
    }

    return f;
}

double Interpolator::sampleFromGrayscaleImage(const ImageRAM* tr, const vec2& position)
{
    vec2 locPix(position[0] - int(position[0]), position[1] - int(position[1]));

    auto f00 = tr->readPixel(size2_t(position[0],     position[1]),     LayerType::Color)[0];
    auto f10 = tr->readPixel(size2_t(position[0] + 1, position[1]),     LayerType::Color)[0];
    auto f01 = tr->readPixel(size2_t(position[0],     position[1] + 1), LayerType::Color)[0];
    auto f11 = tr->readPixel(size2_t(position[0] + 1, position[1] + 1), LayerType::Color)[0];

    return (1 - locPix[1]) * ((1 - locPix[0]) * f00 + locPix[0] * f10) +
                locPix[1]  * ((1 - locPix[0]) * f01 + locPix[0] * f11);
}

}  // namespace inviwo
