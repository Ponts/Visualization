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
#include <labtopo/utils/gradients.h>

/*
a) 
x = -y - 3
y = x - 5

b)
x = -x - 2
y = y - 7

c)
x = (-y - 3)*(-x - 2)
y = (x - 5)*( y - 7)

d)
Partial derivatives
x = -sin(x)
y = cos(y)

*/

namespace inviwo
{

const vec4 Topology::ColorsCP[7] =
    {
        vec4(1, 1, 0, 1),  // Saddle
        vec4(0, 0, 1, 1),  // AttractingNode
        vec4(1, 0, 0, 1),  // RepellingNode
        vec4(0.5, 0, 1, 1),// AttractingFocus
        vec4(1, 0.5, 0, 1),// RepellingFocus
        vec4(0, 1, 0, 1),   // Center
		vec4(1,1,1,1)		// Line

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
		mat2 jacobian = Interpolator::sampleJacobian(vol.get(), critPoints[i]);
		types.push_back(identify(jacobian));
		addVertice(critPoints[i], vertices, types[i]);
		indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size() - 1));
		if (types[i] == TypeCP::Saddle) {
			
			computeSeparatrices(critPoints[i], jacobian, vol.get(), mesh, vertices);
		}
		
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
		vec3(0), vec3(0), ColorsCP[static_cast<int>(type)] });
}

Topology::TypeCP Topology::identify(const mat2& jacobian) {
	auto eigenRes = util::eigenAnalysis(jacobian);
	if (abs(eigenRes.eigenvaluesRe[0]) < 0.01 && abs(eigenRes.eigenvaluesRe[1]) < 0.01)
		return TypeCP::Center;
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
	return TypeCP::Saddle;
}

void Topology::computeSeparatrices(const vec2& pos, const mat2& jacobian, const Volume* vr, std::shared_ptr<BasicMesh> mesh, std::vector<BasicMesh::Vertex>& vertices) {
	auto e = util::eigenAnalysis(jacobian);
	vec2 startPos;
	for (auto v = 0; v < 2; ++v) {
		for (auto i = 0; i < 2; ++i)
			startPos[i] = pos[i] + 0.05*e.eigenvectors[v][i];
		int dir = (e.eigenvaluesRe[v] > 0) ? 1 : -1;
		auto buffer = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
		addVertice(pos, vertices, TypeCP::Line);
		buffer->add(static_cast<std::uint32_t>(vertices.size() - 1));
		integrateSeparatrice(vr, buffer, vertices, startPos, dir);
		for (auto i = 0; i < 2; ++i)
			startPos[i] = pos[i] - 0.05*e.eigenvectors[v][i];
		buffer = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
		addVertice(pos, vertices, TypeCP::Line);
		buffer->add(static_cast<std::uint32_t>(vertices.size() - 1));
		integrateSeparatrice(vr, buffer, vertices, startPos, dir);
	}
	
	

}

void Topology::integrateSeparatrice(const Volume* vr, IndexBufferRAM* buffer, std::vector<BasicMesh::Vertex>& vertices, const vec2& pos, const int dir) {
	addVertice(pos, vertices, TypeCP::Line);
	buffer->add(static_cast<std::uint32_t>(vertices.size() - 1));
	vec2 current = pos;
	double stepsize = 0.05;
	for (int i = 0; i < 200; i++) {
		vec2 previous = current;
		current = Integrator::RK4(vr, current, stepsize, dir);
		double newDist = Integrator::vecLength(current - previous);
		//draw
		addVertice(current, vertices, TypeCP::Line);
		buffer->add(static_cast<std::uint32_t>(vertices.size() - 1));
		if ((current[0] < 0 || current[0] > dims[0] - 1 || current[1] < 0 || current[1] > dims[1] - 1))
			return;
		if (newDist / stepsize < 0.01)
			return;
	}
}


}// namespace
