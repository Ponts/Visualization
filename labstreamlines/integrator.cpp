/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/integrator.h>

namespace inviwo
{

Integrator::Integrator()
{
}

vec2 Integrator::sampleFromField(const VolumeRAM* vr, size3_t dims, const vec2& position)
{
    // Sampled outside the domain!
    if (position[0] < 0 || position[0] > dims[0] - 1 || position[1] < 0 || position[1] > dims[1] - 1)
    {
        return vec2(0, 0);
    }

    int x0 = int(position[0]);
    int y0 = int(position[1]);

    // Leads to accessing only inside the volume
    // Coefficients computation takes care of using the correct values
    if (x0 == dims[0] - 1)
    {
        x0--;
    }
    if (y0 == dims[1] - 1)
    {
        y0--;
    }

    auto f00 = vr->getAsDVec2(size3_t(x0, y0, 0));
    auto f10 = vr->getAsDVec2(size3_t(x0 + 1, y0, 0));
    auto f01 = vr->getAsDVec2(size3_t(x0, y0 + 1, 0));
    auto f11 = vr->getAsDVec2(size3_t(x0 + 1, y0 + 1, 0));

    float x = position[0] - x0;
    float y = position[1] - y0;

    vec2 f;

    for (int i = 0; i < 2; i++)
    {
        f[i] = f00[i] * (1 - x) * (1 - y) + f01[i] * (1 - x) * y + f10[i] * x * (1 - y) + f11[i] * x * y;
    }

    return f;
}


// TODO: Implement a single integration step here

vec2 Integrator::Euler(const VolumeRAM* vr, size3_t dims, const vec2& position, const double stepSize)
{
// Access the vector field with sampleFromField(vr, dims, ...)
	vec2 next;
	vec2 vec = sampleFromField(vr, dims, position);
	for (int i = 0; i < 2; i++) 
		next[i] = position[i] + stepSize * vec[i];

	return next;
}

vec2 Integrator::RK4(const VolumeRAM* vr, size3_t dims, const vec2& position, const double stepSize)
{
	vec2 next;
	vec2 v1 = sampleFromField(vr, dims, position);
	vec2 v2 = sampleFromField(vr, dims, vec2(position.x + 0.5 * stepSize * v1.x, position.y + 0.5 * stepSize * v1.y));
	vec2 v3 = sampleFromField(vr, dims, vec2(position.x + 0.5 * stepSize * v2.x, position.y + 0.5 * stepSize * v2.y));
	vec2 v4 = sampleFromField(vr, dims, vec2(position.x + stepSize * v3.x, position.y + stepSize * v3.y));
	for (int i = 0; i < 2; i++)
		next[i] = position[i] + stepSize*((v1[i] / 6.) + (v2[i] / 3.) + (v3[i] / 3.) + (v4[i] / 6.));
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

vec2 Integrator::RK4(const VolumeRAM* vr, size3_t dims, const vec2& position, const double stepSize, const struct options & opt)
{
	vec2 next;
	vec2 v1 = sampleFromField(vr, dims, position);
	if (opt.normalize)
		normalize(v1);
	vec2 v2 = sampleFromField(vr, dims, vec2(position.x + 0.5 * stepSize * v1.x, position.y + 0.5 * stepSize * v1.y));
	if (opt.normalize)
		normalize(v2);
	vec2 v3 = sampleFromField(vr, dims, vec2(position.x + 0.5 * stepSize * v2.x, position.y + 0.5 * stepSize * v2.y));
	if (opt.normalize)
		normalize(v3);
	vec2 v4 = sampleFromField(vr, dims, vec2(position.x + stepSize * v3.x, position.y + stepSize * v3.y));
	if (opt.normalize)
		normalize(v4);
	for (int i = 0; i < 2; i++)
		next[i] = position[i] + opt.reverse * stepSize*((v1[i] / 6.) + (v2[i] / 3.) + (v3[i] / 3.) + (v4[i] / 6.));
	return next;
}

} // namespace

