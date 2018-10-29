/*********************************************************************
*  Author  : Anke Friederici
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
**********************************************************************/

#pragma once


#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/ports/meshport.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/properties/eventproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <labtopo/labtopomoduledefine.h>



namespace inviwo
{

/** \docpage{org.inviwo.Topology, Vector Field Topology}
    ![](org.inviwo.Topology.png?classIdentifier=org.inviwo.Topology)

    Generate the topological skeleton of a vector field.

    ### Inports
      * __data__ The input here is 2-dimensional vector field (with vectors of
      two components thus two values within each voxel) but it is represented
      by a 3-dimensional volume.
      This processor deals with 2-dimensional data only, therefore it is assumed
      the z-dimension will have size 1 otherwise the 0th slice of the volume
      will be processed.
    

    ### Outports
      * __outMesh__ The output mesh contains points and linesegments that make up
      the topological skeleton.
*/
class IVW_MODULE_LABTOPO_API Topology : public Processor
{
public:
    // All possible first order critical points in 2D vector fields.
    enum class TypeCP
    {
        Saddle = 0,
        AttractingNode = 1,
        RepellingNode = 2,
        AttractingFocus = 3,
        RepellingFocus = 4,
        Center = 5,
		Line = 6,
		Switch = 7
    };

    // Colors according to the TypeCP enum.
    static const vec4 ColorsCP[8];

    // Construction / Deconstruction
  public:
    Topology();
    virtual ~Topology() = default;

    // Methods
  public:
    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;

  protected:
    // Our main computation function
    virtual void process() override;

    //TODO: You may want to declare additional functions here, e.g., extractCriticalPoints.
	bool changeOfSignTest(const Volume* vr, const double x, const double y, const double stepsize);
	bool findCriticalPoint(const Volume* vr, const double x, double y, const double stepsize, std::vector<vec2>& ret);
	void addVertice(const vec2 pos, std::vector<BasicMesh::Vertex>& vertices, const Topology::TypeCP& type);
	Topology::TypeCP identify(const mat2& jacobian);
	void computeSeparatrices(const vec2& pos, const mat2& jacobian, const Volume* vr, std::shared_ptr<BasicMesh> mesh, std::vector<BasicMesh::Vertex>& vertices);
	void integrateSeparatrice(const Volume* vr, IndexBufferRAM* buffer, std::vector<BasicMesh::Vertex>& vertices, const vec2& pos, const int dir);
	void findBoundarySwitchPoints(const Volume* vol, IndexBufferRAM* pointsBuffer, std::vector<BasicMesh::Vertex>& vertices, std::shared_ptr<BasicMesh> mesh);
	int sgn(const double x);
	// Ports
  public:
    // Input data
    VolumeInport inData;

    // Output mesh
    MeshOutport outMesh;
	IntProperty numOfSteps;
	BoolProperty doBoundarySwitch;
  private:
	uvec3 dims;
};

}// namespace inviwo
