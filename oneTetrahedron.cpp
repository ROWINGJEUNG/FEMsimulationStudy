/*sofa framework �̿��� ������� �� ���� ���
1. gmash.exe ���α׷����� msh ������ �����. (mechanicalobject���� dof�� �������� �ʱ� ���� �ʿ���)
2. mechanicalobject���� 4������ �̻��� �����Ѵ�.
3. FFD�� ���� finite cell method�� �����Ѵ�
--> ������ ������� �����ϴ� ���� �ʿ��ϴ�.*/

#include <sofa/core/objectmodel/Context.h>
#include <sofa/core/VecId.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>

#include <sofa/gui/GUIManager.h>
#include <sofa/gui/Main.h>

#include <sofa/helper/ArgumentParser.h>
#include <sofa/helper/UnitTest.h>
#include <sofa/helper/vector_algebra.h>
#include <sofa/helper/vector.h>
#include <sofa/helper/BackTrace.h>
#include <sofa/helper/system/PluginManager.h>
#include <sofa/helper/system/FileRepository.h>

#include <SofaComponentBase/initComponentBase.h>
#include <SofaComponentCommon/initComponentCommon.h>
#include <SofaComponentGeneral/initComponentGeneral.h>
#include <SofaComponentAdvanced/initComponentAdvanced.h>
#include <SofaComponentMisc/initComponentMisc.h>
#include <SofaMiscMapping/SubsetMultiMapping.h>

#include <SofaBaseLinearSolver/CGLinearSolver.h>
#include <SofaBaseMechanics/BarycentricMapping.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaBaseMechanics/UniformMass.h>

#include <SofaBaseTopology/MeshTopology.h>
#include <SofaBaseTopology/EdgeSetTopologyContainer.h>
#include <SofaBaseTopology/RegularGridTopology.h>
#include <SofaBaseTopology/SparseGridTopology.h>

#include <SofaBaseCollision/SphereModel.h>
#include <SofaGeneralTopology/CubeTopology.h>

#include <SofaBaseVisual/VisualStyle.h>
#include <SofaBoundaryCondition/FixedConstraint.h>
#include <SofaImplicitOdeSolver/EulerImplicitSolver.h>
#include <SofaOpenglVisual/OglModel.h>
#include <SofaLoader\MeshObjLoader.h>
#include <SofaGeneralLoader/MeshGmshLoader.h>

#include <SofaSimpleFem/TetrahedronFEMForceField.h>
#include <SofaSimpleFem/HexahedronFEMForceField.h>

#include <SofaDeformable/StiffSpringForceField.h>
#include <SofaDeformable/MeshSpringForceField.h>

#include <sofa/simulation/Node.h>
#include <sofa/simulation/Simulation.h>

#include <SofaSimulationTree/init.h>
#include <SofaSimulationTree/GNode.h>
#include <SofaSimulationTree/TreeSimulation.h>
#include <SofaSimulationGraph/SimpleApi.h>

#include <SceneCreator/SceneCreator.h>

#include <iostream>
#include <sstream>
#include <fstream>

using sofa::core::VecId;
using sofa::defaulttype::Vec3Types;
using sofa::defaulttype::ExtVec3Types;
using Deriv3 = sofa::defaulttype::Vec3Types::Deriv;
using Coord3 = sofa::defaulttype::Vector3;
using VecCoord3 = sofa::helper::vector<Coord3>;


using sofa::core::objectmodel::Data;
using sofa::core::objectmodel::New;
using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;
using sofa::simulation::Node;

using CGLinearSolver = sofa::component::linearsolver::CGLinearSolver<sofa::component::linearsolver::GraphScatteredMatrix, sofa::component::linearsolver::GraphScatteredVector>;
using BarycentricMapping3 = sofa::component::mapping::BarycentricMapping<Vec3Types, Vec3Types>;
using MechanicalObject3 = sofa::component::container::MechanicalObject<Vec3Types>;
using UniformMass3 = sofa::component::mass::UniformMass<Vec3Types, SReal>;


using sofa::component::topology::MeshTopology;
using sofa::component::topology::SparseGridTopology;
using sofa::component::topology::MeshTopology;
using sofa::component::topology::SparseGridTopology;

using sofa::component::visualmodel::VisualStyle;
using FixedConstraint3 = sofa::component::projectiveconstraintset::FixedConstraint<Vec3Types>;
using sofa::component::odesolver::EulerImplicitSolver;
using sofa::component::visualmodel::OglModel;
using sofa::component::loader::MeshObjLoader;
using sofa::component::loader::MeshGmshLoader;

using TetrahedronFEMForceField3 = sofa::component::forcefield::TetrahedronFEMForceField<Vec3Types>;
using HexahedronFEMForceField3 = sofa::component::forcefield::HexahedronFEMForceField<Vec3Types>;
using StiffSpringForceField3 = sofa::component::interactionforcefield::StiffSpringForceField<Vec3Types>;

using sofa::simpleapi::str;
using sofa::simpleapi::createObject;
using sofa::simpleapi::createChild;


int main(int argc, char** argv)
{
	std::vector<std::string> collisionModelTypes = { "Triangle", "Line", "Point" };

    sofa::simulation::tree::init();
	//sofa::helper::BackTrace::autodump();
	//sofa::core::ExecParams::defaultInstance()->setAspectID(0);
    ArgumentParser argParser(argc, argv);
    sofa::gui::GUIManager::RegisterParameters(&argParser);
    argParser.parse();
    sofa::gui::initMain();
    sofa::gui::GUIManager::Init(argv[0]);

	// Gui ������ ���õ� �κ�
    sofa::component::initComponentBase();
    sofa::component::initComponentCommon();
    sofa::component::initComponentGeneral();
    sofa::component::initComponentAdvanced();
    sofa::component::initComponentMisc();

    // The graph root node : gravity already exists in a GNode by default, Node�� VTK�� renderer�� ���� ������ �Ѵ�.
    sofa::simulation::setSimulation(new sofa::simulation::tree::TreeSimulation());
    //sofa::simulation::Node::SPtr groot = sofa::simulation::getSimulation()->createNewGraph("root");
	sofa::helper::system::PluginManager::getInstance().loadPlugin("SofaMiscCollision");               // �浹 ������ ���� �� plugin�� �ҷ��;� �Ѵ�.
	Node::SPtr groot = sofa::modeling::createRootWithCollisionPipeline();                             // Collision Pipeline�� ������ groot ����
    groot->setGravity( Coord3(0, -10, 0) );                                                           // �߷� ���� 0���� �ָ� �ܷ��� �������� �Ʒ��� �þ�� �ʴ´�.
	
 //   // ��ü �ڵ忡�� solver�� �ѹ��� ���ǵǾ�� �Ѵ�.
 //   EulerImplicitSolver::SPtr solver = sofa::core::objectmodel::New<EulerImplicitSolver>();
 //   solver->setName("solver");
 //   solver->f_printLog.setValue(false);           // �� �ڵ��� �ǹ̴�?
	//solver->f_rayleighMass.setValue(0.000);       // �� �ڵ��� �ǹ̴�?
	//solver->f_rayleighStiffness.setValue(0.005);  // �� �ڵ��� �ǹ̴�?
 //   groot->addObject(solver);

 //   CGLinearSolver::SPtr linearSolver = New<CGLinearSolver>();
 //   linearSolver->setName("linearSolver");
	//linearSolver->f_maxIter.setValue(1000);                    // �� �ڵ��� �ǹ̴�? : ��� �ݺ��� �� ȸ ���� �����ϴ� �ڵ�� ����. �ǽð����� ���� ������ �� ���� �� ����.
	//linearSolver->f_tolerance.setValue(1E-3);                  // �� �ڵ��� �ǹ̴�? :
	//linearSolver->f_smallDenominatorThreshold.setValue(1E-3);  // �� �ڵ��� �ǹ̴�? :
 //   groot->addObject(linearSolver);

	Node::SPtr  chain = groot->createChild("Chain");

		// kidney 3D model�� FEM�� �����ϱ� ���� Node
		// Node::SPtr kidneyFEMNode = groot.get()->createChild("kidneyFEMNode");
		Node::SPtr kidneyFEMNode = sofa::modeling::createEulerSolverNode(chain, "kidneyFEMNode");        
		// Tetrahedron degrees of freedom
		/*DOF �ۼ� ���̵�
		1. resize�� 4 ���� ��(0~3�� ���α׷��� �������� ����), �� ������ �յ��� ������ ���� ä�� �����Ǿ� ����.
		2. ��� ���� ������ ª�������� ������ Ŀ����. (�������Ǵ� ������ Ŀ���� �������� ����)
		3. .msh ������ ����ϸ� dof �ϳ��� �۵� ������ �� ����*/
		MechanicalObject3::SPtr DOF = sofa::core::objectmodel::New<MechanicalObject3>();
		kidneyFEMNode->addObject(DOF);
		DOF->resize(4);
		DOF->setName("DOF");
		WriteAccessor<Data<VecCoord3> > x = *DOF->write(VecId::position());
		x[0] = Coord3(0, 10, 0);
		x[1] = Coord3(10, 0, 0);
		x[2] = Coord3(-10 * 0.5, 0, 10 * 0.866);
		x[3] = Coord3(-10 * 0.5, 0, -10 * 0.866);
		DOF->showObject.setValue(true);
		DOF->showObjectScale.setValue(10.);

		// Tetrahedron uniform mass
		UniformMass3::SPtr mass = sofa::core::objectmodel::New<UniformMass3>();
		kidneyFEMNode->addObject(mass);
		mass->setMass(2);
		mass->setName("mass");

		// Tetrahedron constraints
		FixedConstraint3::SPtr constraints = sofa::core::objectmodel::New<FixedConstraint3>();
		kidneyFEMNode->addObject(constraints);
		constraints->setName("constraints");
		constraints->addConstraint(0);


		// Tetrahedron topology
		MeshTopology::SPtr topology = sofa::core::objectmodel::New<MeshTopology>();
		kidneyFEMNode->addObject(topology);
		topology->setName("topology");
		topology->addTetra(0, 1, 2, 3);
		//topology->addHexa(0, 1, 2, 3, 4, 5, 6, 7);

		TetrahedronFEMForceField3::SPtr kidney_fem = sofa::core::objectmodel::New<TetrahedronFEMForceField3>();
		kidneyFEMNode->addObject(kidney_fem);
		kidney_fem->setName("kidney_fem");
		kidney_fem->setMethod("polar");                      // FEM ����� ������� �״�� �θ� �ȴ�.
		kidney_fem->setUpdateStiffnessMatrix(true);
		kidney_fem->setYoungModulus(30);                     // �������� Ư�� ���� (��ȭ���� �� �����ֱ� ���� 6���� ����)
		kidney_fem->setPoissonRatio(0.45);                   // �������� Ư�� ����

			// kidney 3D model�� visualization �ϱ� ���� Node (visual model�� mapping�� addObject ���־�� �Ѵ�.)
			/*visual model�� ���� ���̵�
			1. visual model�� collision model�� ������ �ִ� ��� visual model�� �޽��� �Ĳ��� ������ �ִ� ���̰� collision model�� �޽��� �����ϰ� ������ �ִ� ����.
			   �� ���� ����� �����ϴ�. �������� �����ϰ� ����� ���� visual model�� collision model�� ���� ����ϴ� ������ ����.*/
			Node::SPtr kidney_model_node = kidneyFEMNode.get()->createChild("kidney_model_node");

			// obj file load and set visual model
			MeshObjLoader::SPtr kidney_obj = sofa::core::objectmodel::New<MeshObjLoader>();        // ���߿� Mesh VTK Loader�� �����Ѵ�.
			kidney_obj->setName("kidney_obj");
			kidney_obj->setFilename("model/kidney.obj");
			kidney_obj->load();
			OglModel::SPtr kidney_visual = sofa::core::objectmodel::New<OglModel>();
			kidney_visual->setName("kidney_visual");
			kidney_visual->setSrc("", kidney_obj.get());
			kidney_model_node->addObject(kidney_visual);

			// The mapping between the tetrahedron (DOF) and the liver (visual)
			BarycentricMapping3::SPtr kidney_mapping = sofa::core::objectmodel::New<BarycentricMapping3>();  // mapping ����� �� ����� �������� ����Ѵ�
			kidney_mapping->setModels(DOF.get(), kidney_visual.get());
			kidney_mapping->setName("kidney_mapping");
			kidney_model_node->addObject(kidney_mapping);

		// make collision node for kidney
		/*collision�� ���� ���̵�
		1. �浹�ϴ� �� ����� ������ �� �پ� ������ ������ �߻��� �� �ִ�.*/
		sofa::modeling::createCollisionNodeVec3(kidneyFEMNode, DOF.get(), "model/kidney_collision.obj", collisionModelTypes);


		// mass 3D model�� FEM�� �����ϱ� ���� Node, �ƹ� ����̳� �� ���� ������� �׽�Ʈ �غ���!!
		//Node::SPtr massFEMNode = groot.get()->createChild("massFEMNode");
		Node::SPtr massFEMNode = sofa::modeling::createEulerSolverNode(chain, "massFEMNode");

		MeshGmshLoader::SPtr  loaderFEM = New<MeshGmshLoader>();
		loaderFEM->setFilename("model/test7.msh");
		loaderFEM->load();
		massFEMNode->addObject(loaderFEM);

		MeshTopology::SPtr meshTorusFEM = sofa::core::objectmodel::New<MeshTopology>();
		meshTorusFEM->setSrc("", loaderFEM.get());
		massFEMNode->addObject(meshTorusFEM);

		MechanicalObject3::SPtr dofFEM = sofa::core::objectmodel::New<MechanicalObject3>(); dofFEM->setName("FEM Object");
		massFEMNode->addObject(dofFEM);

		UniformMass3::SPtr uniMassFEM = sofa::core::objectmodel::New<UniformMass3>();
		uniMassFEM->setTotalMass(5); //the whole object will have 5 as given mass
		massFEMNode->addObject(uniMassFEM);

		TetrahedronFEMForceField3::SPtr tetraFEMFF = sofa::core::objectmodel::New<TetrahedronFEMForceField3>();
		tetraFEMFF->setName("FEM");
		tetraFEMFF->setComputeGlobalMatrix(false);
		tetraFEMFF->setMethod("large");
		tetraFEMFF->setPoissonRatio(0.3);
		tetraFEMFF->setYoungModulus(1000);
		massFEMNode->addObject(tetraFEMFF);

		sofa::modeling::createVisualNodeVec3(massFEMNode, dofFEM.get(), "model/mass.obj", "red");

		//// make collision node for tumor -> ���� ���ϴ� ��� �۵��ϴ� �� ������ �ʴ�.
		sofa::modeling::createCollisionNodeVec3(massFEMNode, dofFEM.get(), "model/mass.obj", collisionModelTypes);


		// FEM ���� (���� ���� �׽�Ʈ �غ��� ����, ���� �׽�Ʈ �ʼ�!!)
		{
			//Node::SPtr  torusFEM = sofa::modeling::createEulerSolverNode(chain, "FEM");

			//MeshGmshLoader::SPtr  loaderFEM = New<MeshGmshLoader>();
			//loaderFEM->setFilename(sofa::helper::system::DataRepository.getFile("mesh/torus_low_res.msh"));
			//loaderFEM->load();
			//torusFEM->addObject(loaderFEM);

			//MeshTopology::SPtr meshTorusFEM = sofa::core::objectmodel::New<MeshTopology>();
			//meshTorusFEM->setSrc("", loaderFEM.get());
			//torusFEM->addObject(meshTorusFEM);

			//const Deriv3 translation(2.5, 0, 0);
			//const Deriv3 rotation(90, 0, 0);

			//MechanicalObject3::SPtr dofFEM = sofa::core::objectmodel::New<MechanicalObject3>(); dofFEM->setName("FEM Object");
			//dofFEM->setTranslation(translation[0], translation[1], translation[2]);
			//dofFEM->setRotation(rotation[0], rotation[1], rotation[2]);
			//torusFEM->addObject(dofFEM);

			//UniformMass3::SPtr uniMassFEM = sofa::core::objectmodel::New<UniformMass3>();
			//uniMassFEM->setTotalMass(5); //the whole object will have 5 as given mass
			//torusFEM->addObject(uniMassFEM);

			//TetrahedronFEMForceField3::SPtr tetraFEMFF = sofa::core::objectmodel::New<TetrahedronFEMForceField3>();
			//tetraFEMFF->setName("FEM");
			//tetraFEMFF->setComputeGlobalMatrix(false);
			//tetraFEMFF->setMethod("large");
			//tetraFEMFF->setPoissonRatio(0.3);
			//tetraFEMFF->setYoungModulus(1000);
			//torusFEM->addObject(tetraFEMFF);

			//// Visual node
			//sofa::modeling::createVisualNodeVec3(torusFEM, dofFEM.get(), visualModel, "red", translation, rotation);

			//// Collision node
			//sofa::modeling::createCollisionNodeVec3(torusFEM, dofFEM.get(), collisionModel, collisionModelTypes, translation, rotation);
		}

		// FFD ���� (���� ���� �׽�Ʈ �غ��� ����, ���� �׽�Ʈ �ʼ�!!, FCM�� ������ ������� ����)
		{
			//Node::SPtr  torusFFD = sofa::modeling::createEulerSolverNode(chain, "FFD");

			//const Deriv3 translation(7.5, 0, 0);
			//const Deriv3 rotation(90, 0, 0);

			//MechanicalObject3::SPtr dofFFD = sofa::core::objectmodel::New<MechanicalObject3>(); dofFFD->setName("FFD Object");
			//dofFFD->setTranslation(translation[0], translation[1], translation[2]);
			//dofFFD->setRotation(rotation[0], rotation[1], rotation[2]);
			//torusFFD->addObject(dofFFD);

			//UniformMass3::SPtr uniMassFFD = sofa::core::objectmodel::New<UniformMass3>();
			//uniMassFFD->setTotalMass(5); //the whole object will have 5 as given mass
			//torusFFD->addObject(uniMassFFD);

			//RegularGridTopology::SPtr gridTopo = sofa::core::objectmodel::New<RegularGridTopology>(6, 2, 5); //dimension of the grid
			//gridTopo->setPos(
			//	-2.5, 2.5,  //Xmin, Xmax
			//	-0.5, 0.5,  //Ymin, Ymax
			//	-2, 2       //Zmin, Zmax
			//);
			//torusFFD->addObject(gridTopo);

			//RegularGridSpringForceField3::SPtr FFDFF = sofa::core::objectmodel::New<RegularGridSpringForceField3>();
			//FFDFF->setName("Springs FFD");
			//FFDFF->setStiffness(200);
			//FFDFF->setDamping(0);
			//torusFFD->addObject(FFDFF);

			//// Visual node
			//sofa::modeling::createVisualNodeVec3(torusFFD, dofFFD.get(), visualModel, "yellow");

			//// Collision node
			//sofa::modeling::createCollisionNodeVec3(torusFFD, dofFFD.get(), collisionModel, collisionModelTypes);
		}

	//// ---------------- Interaction force between the kidney and vessel
	//StiffSpringForceField3::SPtr iff = sofa::core::objectmodel::New<StiffSpringForceField3>(DOF.get(), DOF2.get());
	//iff->setPathObject1("@" + kidneyFEMNode->getName() + "/" + DOF->getName());
	//iff->setPathObject2("@" + massFEMNode->getName() + "/" + DOF2->getName());
	//groot->addObject(iff);
	//iff->setName("F13");
	//iff->addSpring(1, 0, 100., 1., 1);  // ������ 1�� splength�� �ǹ��Ѵ� (splength�� �����ΰ�?)


    // Display Flags
    VisualStyle::SPtr style = sofa::core::objectmodel::New<VisualStyle>();
    groot->addObject(style);
    sofa::core::visual::DisplayFlags& flags = *style->displayFlags.beginEdit();
    flags.setShowNormals(false);
    flags.setShowInteractionForceFields(false);
    flags.setShowMechanicalMappings(false);
    flags.setShowCollisionModels(false);
    flags.setShowBoundingCollisionModels(false);
    flags.setShowMappings(false);
    flags.setShowForceFields(false);
    flags.setShowWireFrame(true);
    flags.setShowVisualModels(true);
    flags.setShowBehaviorModels(true);
    style->displayFlags.endEdit();

    // Init the scene
    sofa::simulation::tree::getSimulation()->init(groot.get());
    groot->setAnimate(false);


    //=======================================
    // Run the main loop
    sofa::gui::GUIManager::MainLoop(groot);

    sofa::simulation::tree::cleanup();
    return 0;
}
