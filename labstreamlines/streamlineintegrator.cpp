/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/streamlineintegrator.h>
#include <inviwo/core/util/utilities.h>
#include <inviwo/core/interaction/events/mouseevent.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",  // Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
	: Processor()
	, inData("volIn")
	, outMesh("meshOut")
	, propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(0), vec2(1024), vec2(0.5))
	, propSeedMode("seedMode", "Seeds")
	, rkNIntegration("rkNIntegration", "Step Amount RK", 30, 0, 300)
	, rkStepSize("rkStepsize", "Step size RK", 0.3, 0, 1)
	, reverse("reverse", "reverse integration", false)
	, normalize("normalize", "normalize vector field", false)
	, stopArcLength("arclengthstop", "Stop after arc length", true)
	, arcLength("maxLengthArc", "Max length of the arc", 3, 0, 10)
	, stopAtBoundary("stopAtBoundary", "Stop if outside boundary", true)
	, stopAtZero("stopAtZero", "Stop if velocity is zero", true)
	, stopAtLowVel("stopAtLowVel", "Stop if velocity is low", false)
	, lowVel("lowvel", "Low velocity to stop", 0.0001, 0, 10)
	, showPoints("showpoints", "Show Points", false)
	, nLines("nLines", "Number of lines", 10, 0, 1000)
    // TODO: Initialize additional properties
    // propertyName("propertyIdentifier", "Display Name of the Propery",
    // default value (optional), minimum value (optional), maximum value (optional), increment
    // (optional)); propertyIdentifier cannot have spaces
    , mouseMoveStart("mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
                     MouseButton::Left, MouseState::Press | MouseState::Move) {
    // Register Ports
    addPort(inData);
    addPort(outMesh);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
    addProperty(propSeedMode);
    addProperty(propStartPoint);
    addProperty(mouseMoveStart);

	addProperty(rkStepSize);
	addProperty(rkNIntegration);

	addProperty(reverse);
	addProperty(normalize);
	addProperty(stopArcLength);
	addProperty(arcLength);
	addProperty(stopAtBoundary);
	addProperty(stopAtZero);
	addProperty(stopAtLowVel);
	addProperty(lowVel);
	addProperty(showPoints);
	addProperty(nLines);

    // TODO: Register additional properties
    // addProperty(propertyName);
	util::hide(arcLength);
	util::hide(lowVel);
	util::hide(nLines);
    // You can hide and show properties for a single seed and hide properties for multiple seeds (TODO)
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart);
			util::hide(nLines);
        } else {
            util::hide(propStartPoint, mouseMoveStart);
			util::show(nLines);
        }
    });
	stopArcLength.onChange([this]() {
		if (stopArcLength.get())
			util::show(arcLength);
		else
			util::hide(arcLength);
	});
	

	stopAtLowVel.onChange([this]() {
		if (stopAtLowVel.get())
			util::show(lowVel);
		else
			util::hide(lowVel);
	});

}

void StreamlineIntegrator::eventMoveStart(Event* event) {
    // Handle mouse interaction only if we
    // are in the mode with a single point
    if (propSeedMode.get() == 1) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();
    // Denormalize to volume dimensions
    mousePos.x *= dims.x - 1;
    mousePos.y *= dims.y - 1;
    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void StreamlineIntegrator::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
	srand(29871982312);
	
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    auto vr = vol->getRepresentation<VolumeRAM>();
    dims = vol->getDimensions();
    // The start point should be inside the volume (set maximum to the upper right corner)
    propStartPoint.setMaxValue(vec2(dims.x - 1, dims.y - 1));
	vec2 startPoint = propStartPoint.get();
	
	Integrator::options opts = { reverse.get() ? -1 : 1, normalize.get() };

    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
	
    if (propSeedMode.get() == 0) {
		auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
		auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
		streamLine(vr, indexBufferPoints, indexBufferLines, startPoint, vertices, opts);
        // Draw start point
        
	}
	else {
		// TODO: Seed multiple stream lines either randomly or using a uniform grid
		// (TODO: Bonus, sample randomly according to magnitude of the vector field)

		for (int i = 0; i < nLines.get(); i++) {
			double startX =fRand(0, dims.x - 1);
			double startY = fRand(0, dims.y - 1);
			auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
			auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
			vec2 start = vec2(startX, startY);
			streamLine(vr, indexBufferPoints, indexBufferLines, start, vertices, opts);
		}
		
		

    }

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

double StreamlineIntegrator::fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

void StreamlineIntegrator::streamLine(const VolumeRAM* vr, IndexBufferRAM* bufferPoints, IndexBufferRAM* bufferLines, const vec2 startPoint, std::vector<BasicMesh::Vertex>& vertices, const Integrator::options& opts) {
	
	vertices.push_back({ vec3(startPoint.x / (dims.x - 1), startPoint.y / (dims.y - 1), 0),
		vec3(0), vec3(0), vec4(0, 0, 0, 1) });
	if (showPoints.get())
		bufferPoints->add(static_cast<std::uint32_t>(vertices.size()-1));
	bufferLines->add(static_cast<std::uint32_t>(vertices.size()-1));
	// TODO: Create one stream line from the given start point
	vec2 current = startPoint;
	double distance = 0.;
	for (int i = 0; i < rkNIntegration.get(); i++) {
		vec2 previous = current;
		current = Integrator::RK4(vr, dims, current, rkStepSize.get(), opts);
		double newDist = Integrator::vecLength(current - previous);
		distance += newDist;
		//draw
		vertices.push_back({ vec3(current.x / (dims.x - 1), current.y / (dims.y - 1), 0),
			vec3(0), vec3(0), vec4(0, 0, 255, 1) });
		if (showPoints.get())
			bufferPoints->add(static_cast<std::uint32_t>(vertices.size() - 1));
		bufferLines->add(static_cast<std::uint32_t>(vertices.size() - 1));
		if (stopArcLength.get() && arcLength.get() < distance)
			return;
		if (stopAtBoundary.get() && (current[0] < 0 || current[0] > dims[0] - 1 || current[1] < 0 || current[1] > dims[1] - 1))
			return;
		if (stopAtZero.get() && newDist == 0)
			return;
		if (stopAtLowVel.get() && newDist / rkStepSize.get() < lowVel.get())
			return;



	}
}

}  // namespace inviwo
