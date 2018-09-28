/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <lablic/integrator.h>
#include <lablic/interpolator.h>

namespace inviwo {

Integrator::Integrator() {}

// TODO: Implementation of the functions from last lab
// HINT: There is a change in sampleFromField():
//      Interpolator::sampleFromField(vol.get(), somePosition);
vec2 Integrator::RK4(const Volume* vr, const vec2& position, const double stepSize, const double direction)
{
	vec2 next;
	vec2 v1 = Interpolator::sampleFromField(vr, position);
	normalize(v1);
	vec2 v2 = Interpolator::sampleFromField(vr, vec2(position.x + 0.5 * stepSize * v1.x, position.y + 0.5 * stepSize * v1.y));
	normalize(v2);
	vec2 v3 = Interpolator::sampleFromField(vr, vec2(position.x + 0.5 * stepSize * v2.x, position.y + 0.5 * stepSize * v2.y));
	normalize(v3);
	vec2 v4 = Interpolator::sampleFromField(vr, vec2(position.x + stepSize * v3.x, position.y + stepSize * v3.y));
	normalize(v4);
	for (int i = 0; i < 2; i++)
		next[i] = direction * stepSize*((v1[i] / 6.) + (v2[i] / 3.) + (v3[i] / 3.) + (v4[i] / 6.));
	return next;
}

double Integrator::vecLength(const vec2 & in) {
	return sqrt(in.x*in.x + in.y*in.y);
}

void Integrator::normalize(vec2 & in) {
	double length = sqrt(in.x*in.x + in.y*in.y);
	if (length < 0.001)
		return;
	in.x = in.x / length;
	in.y = in.y / length;
}

}  // namespace inviwo
