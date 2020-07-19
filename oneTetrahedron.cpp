/*sofa framework 이용한 연성장기 모델 구현 방법
1. gmash.exe 프로그램으로 msh 파일을 만든다. (mechanicalobject에서 dof를 정의하지 않기 위해 필요함)
2. mechanicalobject에서 4자유도 이상을 정의한다.
3. FFD를 통해 finite cell method를 구현한다
--> 각각의 장단점을 정리하는 것이 필요하다.*/

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

	// Gui 설정과 관련된 부분
    sofa::component::initComponentBase();
    sofa::component::initComponentCommon();
    sofa::component::initComponentGeneral();
    sofa::component::initComponentAdvanced();
    sofa::component::initComponentMisc();

    // The graph root node : gravity already exists in a GNode by default, Node는 VTK의 renderer와 같은 역할을 한다.
    sofa::simulation::setSimulation(new sofa::simulation::tree::TreeSimulation());
    //sofa::simulation::Node::SPtr groot = sofa::simulation::getSimulation()->createNewGraph("root");
	sofa::helper::system::PluginManager::getInstance().loadPlugin("SofaMiscCollision");               // 충돌 판정을 위해 이 plugin을 불러와야 한다.
	Node::SPtr groot = sofa::modeling::createRootWithCollisionPipeline();                             // Collision Pipeline을 포함한 groot 정의
    groot->setGravity( Coord3(0, -10, 0) );                                                           // 중력 값을 0으로 주면 외력이 없음으로 아래로 늘어나지 않는다.
	
 //   // 전체 코드에서 solver는 한번만 정의되어야 한다.
 //   EulerImplicitSolver::SPtr solver = sofa::core::objectmodel::New<EulerImplicitSolver>();
 //   solver->setName("solver");
 //   solver->f_printLog.setValue(false);           // 이 코드의 의미는?
	//solver->f_rayleighMass.setValue(0.000);       // 이 코드의 의미는?
	//solver->f_rayleighStiffness.setValue(0.005);  // 이 코드의 의미는?
 //   groot->addObject(solver);

 //   CGLinearSolver::SPtr linearSolver = New<CGLinearSolver>();
 //   linearSolver->setName("linearSolver");
	//linearSolver->f_maxIter.setValue(1000);                    // 이 코드의 의미는? : 계산 반복을 몇 회 할지 정의하는 코드로 보임. 실시간성을 위해 제한할 수 있을 것 같음.
	//linearSolver->f_tolerance.setValue(1E-3);                  // 이 코드의 의미는? :
	//linearSolver->f_smallDenominatorThreshold.setValue(1E-3);  // 이 코드의 의미는? :
 //   groot->addObject(linearSolver);

	Node::SPtr  chain = groot->createChild("Chain");

		// kidney 3D model에 FEM을 적용하기 위한 Node
		// Node::SPtr kidneyFEMNode = groot.get()->createChild("kidneyFEMNode");
		Node::SPtr kidneyFEMNode = sofa::modeling::createEulerSolverNode(chain, "kidneyFEMNode");        
		// Tetrahedron degrees of freedom
		/*DOF 작성 가이드
		1. resize는 4 여야 함(0~3은 프로그램이 동작하지 않음), 각 노드들이 균등한 질량을 가진 채로 분포되어 있음.
		2. 노드 사이 간격이 짧아질수록 진폭이 커진다. (단위부피당 질량이 커지기 때문으로 보임)
		3. .msh 파일을 사용하면 dof 하나로 작동 가능할 것 같음*/
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
		kidney_fem->setMethod("polar");                      // FEM 방법과 상관없이 그대로 두면 된다.
		kidney_fem->setUpdateStiffnessMatrix(true);
		kidney_fem->setYoungModulus(30);                     // 물리적인 특성 설정 (변화량을 잘 보여주기 위해 6으로 설정)
		kidney_fem->setPoissonRatio(0.45);                   // 물리적인 특성 설정

			// kidney 3D model을 visualization 하기 위한 Node (visual model과 mapping을 addObject 해주어야 한다.)
			/*visual model에 대한 가이드
			1. visual model과 collision model이 나눠져 있는 경우 visual model은 메쉬가 꼼꼼히 나눠져 있는 모델이고 collision model은 메쉬가 간단하게 나눠져 있는 모델임.
			   두 모델의 모양은 동일하다. 계산과정을 간단하게 만들기 위해 visual model과 collision model로 나눠 사용하는 것으로 보임.*/
			Node::SPtr kidney_model_node = kidneyFEMNode.get()->createChild("kidney_model_node");

			// obj file load and set visual model
			MeshObjLoader::SPtr kidney_obj = sofa::core::objectmodel::New<MeshObjLoader>();        // 나중에 Mesh VTK Loader로 변경한다.
			kidney_obj->setName("kidney_obj");
			kidney_obj->setFilename("model/kidney.obj");
			kidney_obj->load();
			OglModel::SPtr kidney_visual = sofa::core::objectmodel::New<OglModel>();
			kidney_visual->setName("kidney_visual");
			kidney_visual->setSrc("", kidney_obj.get());
			kidney_model_node->addObject(kidney_visual);

			// The mapping between the tetrahedron (DOF) and the liver (visual)
			BarycentricMapping3::SPtr kidney_mapping = sofa::core::objectmodel::New<BarycentricMapping3>();  // mapping 방법은 이 방법을 공통으로 사용한다
			kidney_mapping->setModels(DOF.get(), kidney_visual.get());
			kidney_mapping->setName("kidney_mapping");
			kidney_model_node->addObject(kidney_mapping);

		// make collision node for kidney
		/*collision에 대한 가이드
		1. 충돌하는 두 대상이 완전히 딱 붙어 있으면 오류가 발생할 수 있다.*/
		sofa::modeling::createCollisionNodeVec3(kidneyFEMNode, DOF.get(), "model/kidney_collision.obj", collisionModelTypes);


		// mass 3D model에 FEM을 적용하기 위한 Node, 아무 기능이나 이 모델을 대상으로 테스트 해보자!!
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

		//// make collision node for tumor -> 내가 원하는 대로 작동하는 것 같지는 않다.
		sofa::modeling::createCollisionNodeVec3(massFEMNode, dofFEM.get(), "model/mass.obj", collisionModelTypes);


		// FEM 예시 (아직 구현 테스트 해보지 못함, 구현 테스트 필수!!)
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

		// FFD 예시 (아직 구현 테스트 해보지 못함, 구현 테스트 필수!!, FCM과 유사한 방법으로 보임)
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
	//iff->addSpring(1, 0, 100., 1., 1);  // 마지막 1은 splength를 의미한다 (splength가 무엇인가?)


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
