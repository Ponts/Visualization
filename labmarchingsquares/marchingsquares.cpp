/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp, ...
 *  Init    : Monday, September 11, 2017 - 12:58:42
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/util/utilities.h>


namespace inviwo
{

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo MarchingSquares::processorInfo_
{
    "org.inviwo.MarchingSquares",      // Class identifier
    "Marching Squares",                // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,           // Code state
    Tags::None,                        // Tags
};

const ProcessorInfo MarchingSquares::getProcessorInfo() const
{
    return processorInfo_;
}


MarchingSquares::MarchingSquares()
    :Processor()
    , inData("volumeIn")
    , meshOut("meshOut")
    , propShowGrid("showGrid", "Show Grid")
    , propDeciderType("deciderType", "Decider Type")
    , propMultiple("multiple", "Iso Levels")
    , propIsoValue("isovalue", "Iso Value")
    , propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f),
        vec4(0.0f), vec4(1.0f), vec4(0.1f),
        InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f),
        vec4(0.0f), vec4(1.0f), vec4(0.1f),
        InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &inData)
{
    // Register ports
    addPort(inData);
    addPort(meshOut);
	
    // Register properties
    addProperty(propShowGrid);
    addProperty(propGridColor);
	
    addProperty(propDeciderType);
    propDeciderType.addOption("midpoint", "Mid Point", 0);
    propDeciderType.addOption("asymptotic", "Asymptotic", 1);

    addProperty(propMultiple);
    
    propMultiple.addOption("single", "Single", 0);
    addProperty(propIsoValue);
    addProperty(propIsoColor);

    propMultiple.addOption("multiple", "Multiple", 1);
    addProperty(propNumContours);
    addProperty(propIsoTransferFunc);

    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clearPoints();
    //propIsoTransferFunc.get().addPoint(vec2(0.0f, 1.0f), vec4(0.0f, 1.0f, 1.0f, 1.0f));
	propIsoTransferFunc.get().addPoint(vec2(0.5f, 1.0f), vec4(0.0f, 1.0f, 0.0f, 1.0f));
    //propIsoTransferFunc.get().addPoint(vec2(1.0f, 1.0f), vec4(1.0f, 0.0f, 0.0f, 1.0f));
    propIsoTransferFunc.setCurrentStateAsDefault();

    util::hide(propGridColor, propNumContours, propIsoTransferFunc);
    // Show the grid color property only if grid is actually displayed
    propShowGrid.onChange([this]()
    {
        if (propShowGrid.get())
        {
            util::show(propGridColor);
        }
        else
        {
            util::hide(propGridColor);
        }
    });

    // Show options based on display of one or multiple iso contours
    propMultiple.onChange([this]()
    {
        if (propMultiple.get() == 0)
        {
            util::show(propIsoValue, propIsoColor);
            util::hide(propNumContours, propIsoTransferFunc);
        }
        else
        {
            //util::hide(propIsoValue);
            //util::show(propIsoColor, propNumContours);
            
            // TODO (Bonus): Comment out above if you are using the transfer function
            // and comment in below instead
            util::hide(propIsoValue, propIsoColor);
            util::show(propNumContours, propIsoTransferFunc);
        }
    });

}

void MarchingSquares::process()
{
    if (!inData.hasData()) {
	    return;
    }

    // This results in a shared pointer to a volume
    auto vol = inData.getData();

    // Extract the minimum and maximum value from the input data
    const double minValue = vol->dataMap_.valueRange[0];
    const double maxValue = vol->dataMap_.valueRange[1];

    // Set the range for the isovalue to that minimum and maximum
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);

    // You can print to the Inviwo console with Log-commands:
    LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue << ".");
    // You can also inform about errors and warnings:
    // LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
    // LogProcessorError("I am letting you know about an error"); // Will print error message in red
    // (There is also LogNetwork...() and just Log...(), these display a different source,
    // LogProcessor...() for example displays the name of the processor in the workspace while
    // Log...() displays the identifier of the processor (thus with multiple processors of the
    // same kind you would not know which one the information is coming from
    
    // Retreive data in a form that we can access it
    const VolumeRAM* vr = vol->getRepresentation< VolumeRAM >();
    const size3_t dims = vol->getDimensions();

    // Initialize mesh and vertices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    // Values within the input data are accessed by the function below
    // It's input is the VolumeRAM from above, the dimensions of the volume
    // and the indeces i and j of the position to be accessed where
    // i is in [0, dims.x-1] and j is in [0, dims.y-1]
	
    float valueat00 = getInputValue(vr, dims, 0, 0);
    LogProcessorInfo("Value at (0,0) is: " << valueat00);
    // You can assume that dims.z = 1 and do not need to consider others cases

    // TODO (Bonus) Gaussian filter
    // Our input is const, but you need to compute smoothed data and write it somewhere
    // Create an editable volume like this:
    // Volume volSmoothed(vol->getDimensions(), vol->getDataFormat());
    // auto vrSmoothed = volSmoothed.getEditableRepresentation<VolumeRAM>();
    // Values can be set with
    // vrSmoothed->setFromDouble(vec3(i,j,0), value);
    // getting values works with an editable volume as well
    // getInputValue(vrSmoothed, dims, 0, 0);

    // Grid
	double xCellSize = dims.x - 1;
	double yCellSize = dims.y - 1;
    // Properties are accessed with propertyName.get() 
    if (propShowGrid.get())
    {
        // TODO: Add grid lines of the given color 

        // The function drawLineSegments creates two vertices at the specified positions, 
        // that are placed into the Vertex vector defining our mesh. 
        // An index buffer specifies which of those vertices should be grouped into to make up lines/trianges/quads.
        // Here two vertices make up a line segment.
        auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        // Draw a line segment from v1 to v2 with a color, the coordinates in the final 
        // image range from 0 to 1 for both x and y
        //vec2 v1 = vec2(0.5, 0.5);
        //vec2 v2 = vec2(0.7, 0.7);
		

        //drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid, vertices);
		
		for (double x = 0.; x <= 1.0; x+=1./xCellSize) {
			for (double y = 0.; y <= 1.0; y+=1./yCellSize) {
				bool cont = true;
				if (x == 0. && y < 1.0) {
					drawLineSegment(vec2(x,y), vec2(x,y+1./yCellSize), propGridColor.get(), indexBufferGrid, vertices);
					cont = false;
				}
				if (y == 0. && x < 1.0) {
					drawLineSegment(vec2(x, y), vec2(x+1./xCellSize, y), propGridColor.get(), indexBufferGrid, vertices);
					cont = false;
				}
				if (cont) {
					vec2 v1 = vec2(x - (1. / xCellSize), y);
					vec2 v2 = vec2(x, y - 1. / yCellSize);
					vec2 v0 = vec2(x, y);
					drawLineSegment(v0, v1, propGridColor.get(), indexBufferGrid, vertices);
					drawLineSegment(v0, v2, propGridColor.get(), indexBufferGrid, vertices);
				}		
			}
		}
    }

    // Iso contours

    if (propMultiple.get() == 0)
    {
		drawIsoLine(mesh, vr, vertices, dims, xCellSize, yCellSize, propIsoValue, propIsoColor.get());

    }
    else
    {
        // TODO: Draw the given number (propNumContours) of isolines between 
        // the minimum and maximum value
		double n = propNumContours.get();
		double stepSize = (maxValue - minValue) / (n + 1);
		for (double i = 1; i <= n; i++) {
			double value = minValue + i*stepSize;
			drawIsoLine(mesh, vr, vertices, dims, xCellSize, yCellSize, value, propIsoTransferFunc.get().sample(i/n));
		}
        
        // TODO (Bonus): Use the transfer function property to assign a color
        // The transfer function normalizes the input data and sampling colors
        // from the transfer function assumes normalized input, that means
        // vec4 color = propIsoTransferFunc.get().sample(0.0f);
        // is the color for the minimum value in the data
        // vec4 color = propIsoTransferFunc.get().sample(1.0f);
        // is the color for the maximum value in the data

    }

    // Note: It is possible to add multiple index buffers to the same mesh,
    // thus you could for example add one for the grid lines and one for
    // each isoline
    // Also, consider to write helper functions to avoid code duplication
    // e.g. for the computation of a single iso contour

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

double MarchingSquares::getInputValue(const VolumeRAM* data, const size3_t dims, 
    const size_t i, const size_t j) {
    // Check if the indices are withing the dimensions of the volume
    if (i < dims.x && j < dims.y) {
        return data->getAsDouble(size3_t(i, j, 0));
    } else {
        LogProcessorError(
            "Attempting to access data outside the boundaries of the volume, value is set to 0");
        return 0;
    }
}

void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
                                      IndexBufferRAM* indexBuffer,
                                      std::vector<BasicMesh::Vertex>& vertices) {
    // Add first vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    // A vertex has a position, a normal, a texture coordinate and a color
    // we do not use normal or texture coordinate, but still have to specify them
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    // Add second vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

void MarchingSquares::drawIsoLine(std::shared_ptr<inviwo::BasicMesh>& mesh, const VolumeRAM* vr, std::vector<BasicMesh::Vertex>& vertices, const size3_t& dims, const double xCellSize, const double yCellSize, const double isoValue, const vec4 isoColor) {
	auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
	// TODO: Draw a single isoline at the specified isovalue (propIsoValue) 
	// and color it with the specified color (propIsoColor)
	for (int x = 0; x < dims.x - 1; x++) {
		for (int y = 0; y < dims.y - 1; y++) {
			// see which case
			std::vector <vec2> points;
			const double valueat00 = getInputValue(vr, dims, x, y);
			const double valueat10 = getInputValue(vr, dims, x + 1, y);
			const double valueat01 = getInputValue(vr, dims, x, y + 1);
			const double valueat11 = getInputValue(vr, dims, x + 1, y + 1);
			//00 -> 10
			if (isoValue < valueat00 != isoValue < valueat10) {
				double xcord = (isoValue - (x + 1)*valueat00 + x*valueat10) / (valueat10 - valueat00);
				points.push_back(vec2(xcord / xCellSize, double(y) / yCellSize));
			}
			if (isoValue < valueat10 != isoValue < valueat11) {
				double ycord = (isoValue - (y + 1)*valueat10 + y*valueat11) / (valueat11 - valueat10);
				points.push_back(vec2(double(x + 1) / xCellSize, ycord / yCellSize));
			}
			if (isoValue < valueat01 != isoValue < valueat11) {
				double xcord = (isoValue - (x + 1)*valueat01 + x*valueat11) / (valueat11 - valueat01);
				points.push_back(vec2(xcord / xCellSize, double(y + 1) / yCellSize));
			}
			if (isoValue < valueat00 != isoValue < valueat01) {
				double ycord = (isoValue - (y + 1)*valueat00 + y*valueat01) / (valueat01 - valueat00);
				points.push_back(vec2(double(x) / xCellSize, ycord / yCellSize));
			}

			if (points.size() == 2) {
				drawLineSegment(points[0], points[1], isoColor, indexBufferGrid, vertices);
			}
			else if (points.size() == 4) {
				//midpoint decider
				if (propDeciderType.get() == 0) {
					double mid = (valueat00 + valueat10 + valueat01 + valueat11) / 4;
					if (isoValue < valueat00 != isoValue < mid) { //00 sign is different from mid sign
						drawLineSegment(points[0], points[3], isoColor, indexBufferGrid, vertices);
						drawLineSegment(points[1], points[2], isoColor, indexBufferGrid, vertices);
					}
					else {
						drawLineSegment(points[0], points[1], isoColor, indexBufferGrid, vertices);
						drawLineSegment(points[2], points[3], isoColor, indexBufferGrid, vertices);
					}
				}
				//asymptotic decider
				else if (propDeciderType.get() == 1) {
					// sort by x
					std::sort(points.begin(), points.end(), [](const vec2& a, const vec2& b) {return a.x < b.x; });
					drawLineSegment(points[0], points[1], propIsoColor.get(), indexBufferGrid, vertices);
					drawLineSegment(points[2], points[3], propIsoColor.get(), indexBufferGrid, vertices);
				}
			}
		}
	}
}

} // namespace
