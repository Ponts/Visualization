/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labtopo/integrator.h>
#include <labtopo/interpolator.h>

namespace inviwo {

Integrator::Integrator() {}

// TODO: Implementation of the functions defined in the header file integrator.h
vec2 Integrator::RK4(const Volume* vr, const vec2& position, const double stepSize, const int dir)
{
	vec2 next;
	vec2 v1 = Interpolator::sampleFromField(vr, position);
	vec2 v2 = Interpolator::sampleFromField(vr, vec2(position.x + 0.5 * stepSize * v1.x, position.y + 0.5 * stepSize * v1.y));
	vec2 v3 = Interpolator::sampleFromField(vr, vec2(position.x + 0.5 * stepSize * v2.x, position.y + 0.5 * stepSize * v2.y));
	vec2 v4 = Interpolator::sampleFromField(vr, vec2(position.x + stepSize * v3.x, position.y + stepSize * v3.y));
	for (int i = 0; i < 2; i++)
		next[i] = position[i] + dir * stepSize*((v1[i] / 6.) + (v2[i] / 3.) + (v3[i] / 3.) + (v4[i] / 6.));
	return next;
}

double Integrator::vecLength(const vec2 & in) {
	return sqrt(in.x*in.x + in.y*in.y);
}
}  // namespace inviwo
