/*********************************************************************
*  Author  : Anke Friederici
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
**********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labtopo/integrator.h>
#include <labtopo/interpolator.h>
#include <labtopo/topology.h>


namespace inviwo
{

const vec4 Topology::ColorsCP[6] =
    {
        vec4(1, 1, 0, 1),  // Saddle
        vec4(0, 0, 1, 1),  // AttractingNode
        vec4(1, 0, 0, 1),  // RepellingNode
        vec4(0.5, 0, 1, 1),// AttractingFocus
        vec4(1, 0.5, 0, 1),// RepellingFocus
        vec4(0, 1, 0, 1)   // Center
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",  // Class identifier
    "Vector Field Topology",// Display name
    "KTH Lab",              // Category
    CodeState::Experimental,// Code state
    Tags::None,             // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const
{
    return processorInfo_;
}

Topology::Topology()
	: Processor(), outMesh("meshOut"), inData("inData")
	, numOfSteps("numofsteps", "Nr of steps", 5, 0, 100)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment (optional));
// propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
	

    // TODO: Register additional properties
    // addProperty(propertyName);
	addProperty(numOfSteps);
}

void Topology::process()
{
    // Get input
    if (!inData.hasData())
    {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VolumeRAM* vr = vol->getRepresentation<VolumeRAM>();
    dims = vr->getDimensions();

    // Initialize mesh, vertices and index buffers for the two streamlines and the
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer),
    // or use several index buffers with connectivity type adjacency.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.
    // You can use your previous integration code (copy it over or call it from <lablic/integrator.h>).

    // Looping through all values in the vector field.
	std::vector<vec2> critPoints;
	std::vector<TypeCP> types;
	LogProcessorInfo("dims: " << dims[0] << " " << dims[1]);
	for (auto y = 0; y < dims[1] - 1; ++y) {
		for (auto x = 0; x < dims[0] - 1; ++x) {
			if (changeOfSignTest(vol.get(), x, y, 1.)) {
				findCriticalPoint(vol.get(), x, y, 1., critPoints);
			}
		}
	}
	for (int i = 0; i < critPoints.size(); ++i) {
		//mat2 jacobian = Interpolator::sampleJacobian(vol.get(), critPoints[i]);
		//auto eigenRes = util::eigenAnalysis(jacobian);
		//types.push_back(identify(eigenRes));
		addVertice(critPoints[i], vertices, TypeCP::AttractingFocus);
		indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size() - 1));
	}
	LogProcessorInfo(critPoints.size());
    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}


bool Topology::findCriticalPoint(const Volume* vr, const double x, const double y, const double stepsize, std::vector<vec2>& ret) {
	bool cont = true;
	double newStepsize = stepsize * 0.5;
	if (stepsize < pow(0.5, numOfSteps)) {
		ret.push_back(vec2(x + newStepsize, y + newStepsize));
		return false;
	}
	if (changeOfSignTest(vr, x, y, newStepsize)) 
		cont = findCriticalPoint(vr, x, y, newStepsize, ret);
	if (cont && changeOfSignTest(vr, x + newStepsize, y, newStepsize)) 
		cont = findCriticalPoint(vr, x + newStepsize, y, newStepsize, ret);
	if (cont && changeOfSignTest(vr, x, y + newStepsize, newStepsize)) 
		cont = findCriticalPoint(vr, x, y + newStepsize, newStepsize, ret);
	if (cont && changeOfSignTest(vr, x + newStepsize, y + newStepsize, newStepsize)) 
		cont = findCriticalPoint(vr, x + newStepsize, y + newStepsize, newStepsize, ret);

	return cont;
} 


bool Topology::changeOfSignTest(const Volume* vr, const double x, const double y, const double stepsize) {
	dvec2 v00 = Interpolator::sampleFromField(vr, vec2(x, y));
	dvec2 v10 = Interpolator::sampleFromField(vr, vec2(x + stepsize, y));
	dvec2 v01 = Interpolator::sampleFromField(vr, vec2(x, y + stepsize));
	dvec2 v11 = Interpolator::sampleFromField(vr, vec2(x + stepsize, y + stepsize));
	dvec2 points[4] = { v00, v10, v01, v11 };
	int xPlusses = 0;
	int yPlusses = 0;
	for (int i = 0; i < 4; ++i) {
		if (0 < points[i].x)
			xPlusses++;
		if (0 < points[i].y)
			yPlusses++;
	}
	if ((xPlusses == 0) || (xPlusses == 4) || (yPlusses == 0) || (yPlusses == 4))
		return false;
	return true;
}

void Topology::addVertice(const vec2 pos, std::vector<BasicMesh::Vertex>& vertices, const Topology::TypeCP& type) {
	// check what type the vertice is

	vertices.push_back({ vec3(pos.x / (dims[0] - 1), pos.y / (dims[1] - 1), 0),
		vec3(0), vec3(0), ColorsCP[0] });

}

Topology::TypeCP Topology::identify(const util::EigenResult& eigenRes) {
	if (eigenRes.eigenvaluesRe[0] > 0 && eigenRes.eigenvaluesRe[1] > 0) {
		if (eigenRes.eigenvaluesIm[0] == 0)
			return TypeCP::RepellingNode;
		return TypeCP::RepellingFocus;
	}
	if (eigenRes.eigenvaluesRe[0] < 0 && eigenRes.eigenvaluesRe[1] < 0) {
		if (eigenRes.eigenvaluesIm[0] == 0) 
			return TypeCP::AttractingNode;
		return TypeCP::AttractingFocus;
	}
	if (eigenRes.eigenvaluesRe[0] == 0 && eigenRes.eigenvaluesRe[1] == 0) 
		return TypeCP::Center;
	return TypeCP::Saddle;
}

}// namespace
