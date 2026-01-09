//
// Filename     : baseReaction.cc
// Description  : A base class describing variable updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2003
// Revision     : $Id: baseReaction.cc,v 1.25 2006/03/18 00:05:14 henrik Exp $
//
#include "baseReaction.h"

#include <vector>

#include "adhocReaction.h"
#include "bending.h"
#include "boolean.h"
#include "calculate.h"
#include "cellTime.h"
#include "centerTriangulation.h"
#include "creation.h"
#include "degradation.h"
#include "dilution.h"
#include "directionReaction.h"
#include "fiberModel.h"
#include "force.h"
#include "grn.h"
#include "growth.h"
#include "growthForce.h"
#include "hypocotyl3D.h"
#include "initiation.h"
#include "massAction.h"
#include "mechanical.h"
#include "mechanicalSpring.h"
#include "mechanicalTRBS.h"
#include "membraneCycling.h"
#include "membraneCyclingAll.h"
#include "network.h"
#include "pressure2D.h"
#include "sisterVertex.h"
#include "transport.h"
#include "turgorGrowth.h"

BaseReaction::~BaseReaction() {}

BaseReaction *
BaseReaction::createReaction(std::vector<double> &paraValue,
                             std::vector<std::vector<size_t> > &indValue,
                             std::string idValue) {
    // Growth related updates
    // namespace WallGrowth (growth.h,growth.cc)
    if (idValue == "WallGrowth::Constant")
        return new WallGrowth::Constant(paraValue, indValue);
    else if (idValue == "WallGrowth::Stress" || idValue == "WallGrowthStress")
        return new WallGrowth::Stress(paraValue, indValue);
    else if (idValue == "WallGrowth::Strain" || idValue == "WallGrowthStrain")
        return new WallGrowth::Strain(paraValue, indValue);
    else if (idValue == "WallGrowth::CenterTriangulation::Constant" ||
             idValue == "CenterTriangulation::WallGrowth::Constant")
        return new WallGrowth::CenterTriangulation::Constant(paraValue, indValue);
    else if (idValue == "WallGrowth::CenterTriangulation::Stress" ||
             idValue == "CenterTriangulation::WallGrowth::Stress" ||
             idValue == "WallGrowthStresscenterTriangulation")
        return new WallGrowth::CenterTriangulation::Stress(paraValue, indValue);
    else if (idValue == "WallGrowth::CenterTriangulation::StressConcentrationHill" ||
             idValue == "CenterTriangulation::WallGrowth::StressConcentrationHill")
        return new WallGrowth::CenterTriangulation::StressConcentrationHill(paraValue, indValue);
    else if (idValue == "WallGrowth::CenterTriangulation::StrainTRBS" ||
             idValue == "CenterTriangulation::WallGrowth::StrainTRBS")
        return new WallGrowth::CenterTriangulation::StrainTRBS(paraValue, indValue);
    else if (idValue == "WallGrowth::CenterTriangulation::StrainTRBSConcentrationHill" ||
             idValue == "CenterTriangulation::WallGrowth::StrainTRBSConcentrationHill")
        return new WallGrowth::CenterTriangulation::StrainTRBSConcentrationHill(paraValue, indValue);
    else if (idValue == "WallGrowth::CenterTriangulation::VectorTRBS" ||
             idValue == "CenterTriangulation::WallGrowth::VectorTRBS")
        return new WallGrowth::CenterTriangulation::VectorTRBS(paraValue, indValue);
    else if (idValue == "WallGrowth::StressSpatial" || idValue == "WallGrowthStressSpatial")
        return new WallGrowth::StressSpatial(paraValue, indValue);
    else if (idValue == "WallGrowth::StressSpatialSingle" || idValue == "WallGrowthStressSpatialSingle")
        return new WallGrowth::StressSpatialSingle(paraValue, indValue);
    else if (idValue == "WallGrowth::StressConcentrationHill" ||
             idValue == "WallGrowthStressConcentrationHill")
        return new WallGrowth::StressConcentrationHill(paraValue, indValue);
    else if (idValue == "WallGrowth::ConstantStressEpidermalAsymmetric" ||
             idValue == "WallGrowthConstantStressEpidermalAsymmetric")
        return new WallGrowth::ConstantStressEpidermalAsymmetric(paraValue, indValue);
    else if (idValue == "WallGrowth::Force")
        return new WallGrowth::Force(paraValue, indValue);

    // namespace TurgorGrowth (turgorGrowth.h(.cc))
    else if (idValue == "TurgorGrowth::WaterVolume" ||
             idValue == "WaterVolumeFromTurgor")
        return new TurgorGrowth::WaterVolume(paraValue, indValue);

    // namespace Dilution (dilution.h(.cc))
    else if (idValue == "Dilution::FromVertexDerivs" ||
             idValue == "DilutionFromVertexDerivs")
        return new Dilution::FromVertexDerivs(paraValue, indValue);

    // Mechanical interactions between vertices
    // mechanicalSpring.h,mechanicalSpring.cc
    else if (idValue == "WallMechanics::Spring" ||
             idValue == "VertexFromWallSpring")
        return new WallMechanics::Spring(paraValue, indValue);
    else if (idValue == "VertexFromWallSpringMTnew")
        return new VertexFromWallSpringMTnew(paraValue, indValue);
    else if (idValue == "VertexFromWallBoundarySpring")
        return new VertexFromWallBoundarySpring(paraValue, indValue);
    else if (idValue == "CenterTriangulation::EdgeSpring" ||
             idValue == "CenterTriangulation::Spring")
        return new CenterTriangulation::EdgeSpring(paraValue, indValue);
    else if (idValue == "VertexFromDoubleWallSpring")
        return new VertexFromDoubleWallSpring(paraValue, indValue);
    else if (idValue == "VertexFromWallSpringSpatial")
        return new VertexFromWallSpringSpatial(paraValue, indValue);
    else if (idValue == "WallMechanics::SpringConcentrationHill" ||
             idValue == "VertexFromWallSpringConcentrationHill")
        return new WallMechanics::SpringConcentrationHill(paraValue, indValue);
    else if (idValue == "WallMechanics::SpringInternalExternalThreshold" ||
             idValue == "SpringInternalExternalThreshold")
        return new WallMechanics::SpringInternalExternalThreshold(paraValue, indValue);
    else if (idValue == "VertexFromWallSpringMT")
        return new VertexFromWallSpringMT(paraValue, indValue);
    else if (idValue == "VertexFromWallSpringMTSpatial")
        return new VertexFromWallSpringMTSpatial(paraValue, indValue);
    else if (idValue == "VertexFromWallSpringMTHistory")
        return new VertexFromWallSpringMTHistory(paraValue, indValue);
    else if (idValue == "WallMechanics::SpringEpidermal" ||
             idValue == "VertexFromEpidermalWallSpring")
        return new WallMechanics::SpringEpidermal(paraValue, indValue);
    else if (idValue == "WallMechanics::SpringEpidermalCell" ||
             idValue == "VertexFromEpidermalCellWallSpring")
        return new WallMechanics::SpringEpidermalCell(paraValue, indValue);
    else if (idValue == "WallMechanics::ViscoElastic")
        return new WallMechanics::ViscoElastic(paraValue, indValue);
    else if (idValue == "VertexFromWallSpringMTConcentrationHill")
        return new VertexFromWallSpringMTConcentrationHill(paraValue, indValue);
    else if (idValue == "VertexFromDoubleWallSpringMTConcentrationHill")
        return new VertexFromDoubleWallSpringMTConcentrationHill(paraValue, indValue);
    else if (idValue == "VertexFromExternalSpring")
        return new VertexFromExternalSpring(paraValue, indValue);
    else if (idValue == "VertexFromExternalSpringFromPerpVertex")
        return new VertexFromExternalSpringFromPerpVertex(paraValue, indValue);
    else if (idValue == "VertexFromExternalSpringFromPerpVertexDynamic")
        return new VertexFromExternalSpringFromPerpVertexDynamic(paraValue, indValue);
    else if (idValue == "cellcellRepulsion")
        return new cellcellRepulsion(paraValue, indValue);
    else if (idValue == "vertexFromSubstrate")
        return new vertexFromSubstrate(paraValue, indValue);

    // Namespace Pressure2D for forces generated perpendicular to edges (along faces)
    // pressure2D.h, pressure2D.cc
    else if (idValue == "Pressure2D::AreaPotential")
        return new Pressure2D::AreaPotential(paraValue, indValue);
    else if (idValue == "Pressure2D::AreaPotentialTri")
        return new Pressure2D::AreaPotentialTri(paraValue, indValue);
    else if (idValue == "Pressure2D::AreaPotentialTriSpatialThreshold")
        return new Pressure2D::AreaPotentialTriSpatialThreshold(paraValue, indValue);
    else if (idValue == "Pressure2D::AreaPotentialTargetArea")
        return new Pressure2D::AreaPotentialTargetArea(paraValue, indValue);

    // Mechanical interactions between vertices
    // mechanical.h,mechanical.cc
    //  Pressure forces (from 2D cells) implemented assuming a CenterTriangulation
    else if (idValue == "CenterTriangulation::VertexFromCellPressure" ||
             idValue == "VertexFromCellPressurecenterTriangulation")
        return new CenterTriangulation::VertexFromCellPressure(paraValue, indValue);
    else if (idValue == "CenterTriangulation::VertexFromCellPressureLinear" ||
             idValue == "VertexFromCellPressurecenterTriangulationLinear")
        return new CenterTriangulation::VertexFromCellPressureLinear(paraValue, indValue);
    else if (idValue == "PerpendicularWallPressure")
        return new PerpendicularWallPressure(paraValue, indValue);

    else if (idValue == "TargetAreaFromPressure")
        return new TargetAreaFromPressure(paraValue, indValue);
    else if (idValue == "VertexFromCellPowerdiagram")
        return new VertexFromCellPowerdiagram(paraValue, indValue);

    // Namespace Pressure3D for forces generated perpendicular to faces
    else if (idValue == "Pressure3D::Constant" ||
             idValue == "VertexFromCellPlane")
        return new Pressure3D::Constant(paraValue, indValue);
    else if (idValue == "Pressure3D::Linear" ||
             idValue == "VertexFromCellPlaneLinear")
        return new Pressure3D::Linear(paraValue, indValue);
    else if (idValue == "Pressure3D::Spatial" ||
             idValue == "VertexFromCellPlaneSpatial")
        return new Pressure3D::Spatial(paraValue, indValue);
    else if (idValue == "Pressure3D::ConcentrationHill" ||
             idValue == "VertexFromCellPlaneConcentrationHill")
        return new Pressure3D::ConcentrationHill(paraValue, indValue);
    else if (idValue == "Pressure3D::Normalized" ||
             idValue == "VertexFromCellPlaneNormalized")
        return new Pressure3D::Normalized(paraValue, indValue);
    else if (idValue == "Pressure3D::NormalizedSpatial" ||
             idValue == "VertexFromCellPlaneNormalizedSpatial")
        return new Pressure3D::NormalizedSpatial(paraValue, indValue);
    else if (idValue == "Pressure3D::SphereCylinder" ||
             idValue == "VertexFromCellPlaneSphereCylinder")
        return new Pressure3D::SphereCylinder(paraValue, indValue);
    else if (idValue == "Pressure3D::SphereCylinderConcentrationHill" ||
             idValue == "VertexFromCellPlaneSphereCylinderConcentrationHill")
        return new Pressure3D::SphereCylinderConcentrationHill(paraValue, indValue);
    else if (idValue == "Pressure3D::Triangular" ||
             idValue == "VertexFromCellPlaneTriangular")
        return new Pressure3D::Triangular(paraValue, indValue);
    else if (idValue == "Pressure3D::CenterTriangulation::Linear" ||
             idValue == "CenterTriangulation::Pressure3D::Linear" ||
             idValue == "VertexFromCellPlaneLinearCenterTriangulation")
        return new Pressure3D::CenterTriangulation::Linear(paraValue, indValue);

    // Forces acting on vertices, collected in namespace Force
    // force.h, force.cc
    // HJ: some yet needs to be moved from mechanical.h
    else if (idValue == "Force::ForceFromPlane")
        return new Force::ForceFromPlane(paraValue, indValue);
    else if (idValue == "Force::Cylinder" ||
             idValue == "CylinderForce")
        return new Force::Cylinder(paraValue, indValue);
    else if (idValue == "Force::SphereCylinder" ||
             idValue == "SphereCylinderForce")
        return new Force::SphereCylinder(paraValue, indValue);
    else if (idValue == "Force::SphereCylinderRadius" ||
             idValue == "SphereCylinderForceFromRadius")
        return new Force::SphereCylinderRadius(paraValue, indValue);
    else if (idValue == "Force::InfiniteWall" ||
             idValue == "InfiniteWallForce")
        return new Force::InfiniteWall(paraValue, indValue);
    else if (idValue == "Force::EpidermalCoordinate")
        return new Force::EpidermalCoordinate(paraValue, indValue);
    else if (idValue == "Force::EpidermalRadial" ||
             idValue == "EpidermalRadialForce")
        return new Force::EpidermalRadial(paraValue, indValue);
    else if (idValue == "Force::IndexRadial" ||
             idValue == "VertexForceOrigoFromIndex")
        return new Force::IndexRadial(paraValue, indValue);
    else if (idValue == "Force::CellIndexRadial" ||
             idValue == "CellForceOrigoFromIndex")
        return new Force::CellIndexRadial(paraValue, indValue);
    else if (idValue == "Force::Axial")
        return new Force::Axial(paraValue, indValue);
    else if (idValue == "Force::Vector" ||
             idValue == "VertexFromForce")
        return new Force::Vector(paraValue, indValue);
    else if (idValue == "Force::VectorLinear" ||
             idValue == "VertexFromForceLinear")
        return new Force::VectorLinear(paraValue, indValue);
    else if (idValue == "Force::Ball" ||
             idValue == "VertexFromBall")
        return new Force::Ball(paraValue, indValue);
    else if (idValue == "Force::Parabolid" ||
             idValue == "VertexFromParabolid")
        return new Force::Parabolid(paraValue, indValue);
    else if (idValue == "Force::ExternalWall" ||
             idValue == "VertexFromExternalWall")
        return new Force::ExternalWall(paraValue, indValue);

    // Growth by adding forces to vertices
    // GrowthForce.h(.cc)
    else if (idValue == "GrowthForce::Radial" ||
             idValue == "MoveVertexRadially")
        return new GrowthForce::Radial(paraValue, indValue);
    else if (idValue == "GrowthForce::CenterTriangulation::Radial" ||
             idValue == "CenterTriangulation::GrowthForce::Radial" ||
             idValue == "MoveVertexRadiallycenterTriangulation")
        return new GrowthForce::CenterTriangulation::Radial(paraValue, indValue);
    else if (idValue == "GrowthForce::centerTriangulation::ForceToCell" ||
             idValue == "centerTriangulation::GrowthForce::ForceToCell")
        return new GrowthForce::CenterTriangulation::ForceToCell(paraValue, indValue);
    else if (idValue == "GrowthForce::EpidermalRadial" ||
             idValue == "MoveEpidermalVertexRadially")
        return new GrowthForce::EpidermalRadial(paraValue, indValue);
    else if (idValue == "GrowthForce::X" ||
             idValue == "MoveVerteX" ||
             idValue == "MoveVertexX")
        return new GrowthForce::X(paraValue, indValue);
    else if (idValue == "GrowthForce::Y" ||
             idValue == "MoveVertexY")
        return new GrowthForce::Y(paraValue, indValue);
    else if (idValue == "GrowthForce::SphereCylinder" ||
             idValue == "MoveVertexSphereCylinder")
        return new GrowthForce::SphereCylinder(paraValue, indValue);

    // centerTriangulation.h (.cc)
    // Reactions related to a center triangulation of cells
    else if (idValue == "CenterTriangulation::Initiate")
        return new CenterTriangulation::Initiate(paraValue, indValue);

    // mechanicalTRBS.h (.cc)
    // Mechanical updates related to triangular (biquadratic) springs
    else if (idValue == "VertexFromTRBS")
        return new VertexFromTRBS(paraValue, indValue);
    else if (idValue == "VertexFromTRBScenterTriangulation")
        return new VertexFromTRBScenterTriangulation(paraValue, indValue);
    else if (idValue == "VertexFromTRBScenterTriangulationConcentrationHill")
        return new VertexFromTRBScenterTriangulationConcentrationHill(paraValue, indValue);
    else if (idValue == "VertexFromTRBSMT")
        return new VertexFromTRBSMT(paraValue, indValue);
    else if (idValue == "VertexFromTRBScenterTriangulationMT")
        return new VertexFromTRBScenterTriangulationMT(paraValue, indValue);
    else if (idValue == "VertexFromTRLScenterTriangulationMT")
        return new VertexFromTRLScenterTriangulationMT(paraValue, indValue);
    else if (idValue == "VertexFromTRBScenterTriangulationConcentrationHillMT")
        return new VertexFromTRBScenterTriangulationConcentrationHillMT(paraValue, indValue);
    else if (idValue == "VertexFromTRBScenterTriangulationMTOpt")
        return new VertexFromTRBScenterTriangulationMTOpt(paraValue, indValue);

    // FiberModel namespace (fiberModel.h(.cc))
    else if (idValue == "FiberModel::General" ||
             idValue == "FiberModel")
        return new FiberModel::General(paraValue, indValue);
    else if (idValue == "FiberModel::Deposition" ||
             idValue == "FiberDeposition")
        return new FiberModel::Deposition(paraValue, indValue);
    // else if (idValue == "FiberModel::Linear")
    // return new FiberModel::Linear(paraValue, indValue);
    // else if (idValue == "FiberModel::Hill")
    // return new FiberModel::Hill(paraValue, indValue);

    // bending.h (.cc)
    else if (idValue == "Bending::NeighborCenter")
        return new Bending::NeighborCenter(paraValue, indValue);
    else if (idValue == "Bending::Angle")
        return new Bending::Angle(paraValue, indValue);
    else if (idValue == "Bending::AngleInitiate")
        return new Bending::AngleInitiate(paraValue, indValue);
    else if (idValue == "Bending::AngleRelax")
        return new Bending::AngleRelax(paraValue, indValue);

    // creation.h,creation.cc
    else if (idValue == "Creation::Zero" ||
             idValue == "CreationZero")
        return new Creation::Zero(paraValue, indValue);
    else if (idValue == "Creation::One" ||
             idValue == "CreationOne")
        return new Creation::One(paraValue, indValue);
    else if (idValue == "Creation::Two" ||
             idValue == "CreationTwo")
        return new Creation::Two(paraValue, indValue);
    else if (idValue == "Creation::Three" ||
             idValue == "CreationThree")
        return new Creation::Three(paraValue, indValue);
    else if (idValue == "Creation::SpatialSphere" ||
             idValue == "CreationSpatialSphere")
        return new Creation::SpatialSphere(paraValue, indValue);
    else if (idValue == "Creation::SpatialCylinder" ||
             idValue == "CreationSpatialCylinder")
        return new Creation::SpatialCylinder(paraValue, indValue);
    else if (idValue == "Creation::SpatialRing" ||
             idValue == "CreationSpatialRing")
        return new Creation::SpatialRing(paraValue, indValue);
    else if (idValue == "Creation::SpatialCoordinate" ||
             idValue == "CreationSpatialCoordinate")
        return new Creation::SpatialCoordinate(paraValue, indValue);
    else if (idValue == "Creation::SpatialPlane" ||
             idValue == "CreationSpatialPlane")
        return new Creation::SpatialPlane(paraValue, indValue);
    else if (idValue == "Creation::FromList" ||
             idValue == "CreationFromList")
        return new Creation::FromList(paraValue, indValue);
    else if (idValue == "Creation::OneGeometric" ||
             idValue == "CreationOneGeometric")
        return new Creation::OneGeometric(paraValue, indValue);
    else if (idValue == "Creation::Sinus" ||
             idValue == "CreationSinus" ||
             idValue == "creationSinus")  // just because it was accepterd...
        return new Creation::Sinus(paraValue, indValue);

    // degradation.h,degradation.cc
    else if (idValue == "Degradation::One" ||
             idValue == "DegradationOne")
        return new Degradation::One(paraValue, indValue);
    else if (idValue == "Degradation::Two" ||
             idValue == "DegradationTwo")
        return new Degradation::Two(paraValue, indValue);
    else if (idValue == "Degradation::N" ||
             idValue == "DegradationN")
        return new Degradation::N(paraValue, indValue);
    else if (idValue == "Degradation::TwoGeometric" ||
             idValue == "DegradationTwoGeometric")
        return new Degradation::TwoGeometric(paraValue, indValue);
    else if (idValue == "Degradation::Hill" ||
             idValue == "DegradationHill")
        return new Degradation::Hill(paraValue, indValue);
    else if (idValue == "Degradation::HillN" ||
             idValue == "DegradationHillN")
        return new Degradation::HillN(paraValue, indValue);
    else if (idValue == "Degradation::OneWall" ||
             idValue == "DegradationOneWall")
        return new Degradation::OneWall(paraValue, indValue);
    else if (idValue == "Degradation::OneBoundary" ||
             idValue == "DegradationOneBoundary")
        return new Degradation::OneBoundary(paraValue, indValue);
    else if (idValue == "Degradation::OneFromList" ||
             idValue == "DegradationOneFromList")
        return new Degradation::OneFromList(paraValue, indValue);

    // grn.h,grn.cc
    else if (idValue == "Grn::Hill" ||
             idValue == "Hill")
        return new Hill(paraValue, indValue);
    else if (idValue == "HillGeneralOne")
        return new HillGeneralOne(paraValue, indValue);
    else if (idValue == "HillGeneralOne_TwoInputs")
        return new HillGeneralOne_TwoInputs(paraValue, indValue);
    else if (idValue == "HillGeneralTwo")
        return new HillGeneralTwo(paraValue, indValue);
    else if (idValue == "HillGeneralThree")
        return new HillGeneralThree(paraValue, indValue);
    else if (idValue == "Grn")
        return new Grn(paraValue, indValue);
    else if (idValue == "Gsrn2")
        return new Gsrn2(paraValue, indValue);

    // transport.h,transport.cc
    // Namespace Diffusion
    else if (idValue == "Diffusion::MembraneSimple" || idValue == "MembraneDiffusionSimple")
        return new MembraneDiffusionSimple(paraValue, indValue);
    else if (idValue == "Diffusion::Simple" || idValue == "DiffusionSimple")
        return new DiffusionSimple(paraValue, indValue);
    else if (idValue == "Diffusion::SimpleOne" || idValue == "DiffusionSimpleOne")
        return new DiffusionSimpleOne(paraValue, indValue);
    else if (idValue == "Diffusion::ConductiveSimple" || idValue == "DiffusionConductiveSimple")
        return new DiffusionConductiveSimple(paraValue, indValue);
    else if (idValue == "Diffusion::2D" || idValue == "Diffusion2D"
	     || idValue == "Diffusion::2d" || idValue == "Diffusion2d" )
        return new Diffusion2d(paraValue, indValue);
    // Namespace Transport
    else if (idValue == "ActiveTransportCellEfflux")
        return new ActiveTransportCellEfflux(paraValue, indValue);
    else if (idValue == "DiffusionActiveTransportCell")
        return new DiffusionActiveTransportCell(paraValue, indValue);
    else if (idValue == "ActiveTransportCellEffluxMM")
        return new ActiveTransportCellEffluxMM(paraValue, indValue);
    else if (idValue == "ActiveTransportWall")
        return new ActiveTransportWall(paraValue, indValue);
    else if (idValue == "InfluxActiveTransportCell")
        return new InfluxActiveTransportCell(paraValue, indValue);

    // network.h,network.cc
    else if (idValue == "LinPolarizationFast")
        return new LinPolarizationFast(paraValue, indValue);
    else if (idValue == "LinPolarizationFastExact")
        return new LinPolarizationFastExact(paraValue, indValue);
    else if (idValue == "CellCellAuxinTransport")
        return new CellCellAuxinTransport(paraValue, indValue);
    else if (idValue == "AuxinModelSimple1")
        return new AuxinModelSimple1(paraValue, indValue);
    else if (idValue == "AuxinModel1")
        return new AuxinModel1(paraValue, indValue);
    else if (idValue == "AuxinModel1S")
        return new AuxinModel1S(paraValue, indValue);
    else if (idValue == "AuxinModelStress")
        return new AuxinModelStress(paraValue, indValue);
    else if (idValue == "AuxinModelSimpleStress")
        return new AuxinModelSimpleStress(paraValue, indValue);
    else if (idValue == "AuxinModelSimple1Wall")
        return new AuxinModelSimple1Wall(paraValue, indValue);
    else if (idValue == "AuxinModelSimple2")
        return new AuxinModelSimple2(paraValue, indValue);
    else if (idValue == "AuxinModelSimple3")
        return new AuxinModelSimple3(paraValue, indValue);
    else if (idValue == "AuxinModelSimple4")
        return new AuxinModelSimple4(paraValue, indValue);
    else if (idValue == "AuxinModelSimple5")
        return new AuxinModelSimple5(paraValue, indValue);
    else if (idValue == "AuxinModel4")
        return new AuxinModel4(paraValue, indValue);
    else if (idValue == "AuxinModel5")
        return new AuxinModel5(paraValue, indValue);
    else if (idValue == "AuxinModel6")
        return new AuxinModel6(paraValue, indValue);
    else if (idValue == "AuxinModel7")
        return new AuxinModel7(paraValue, indValue);
    else if (idValue == "AuxinTransportCellCellNoGeometry")
        return new AuxinTransportCellCellNoGeometry(paraValue, indValue);
    else if (idValue == "AuxinWallModel")
        return new AuxinWallModel(paraValue, indValue);
    else if (idValue == "AuxinROPModel")
        return new AuxinROPModel(paraValue, indValue);
    else if (idValue == "AuxinROPModel2")
        return new AuxinROPModel2(paraValue, indValue);
    else if (idValue == "AuxinROPModel3")
        return new AuxinROPModel3(paraValue, indValue);
    else if (idValue == "AuxinPINBistabilityModel")
        return new AuxinPINBistabilityModel(paraValue, indValue);
    else if (idValue == "AuxinPINBistabilityModelCell")
        return new AuxinPINBistabilityModelCell(paraValue, indValue);
    else if (idValue == "AuxinExoBistability")
        return new AuxinExoBistability(paraValue, indValue);
    else if (idValue == "AuxinPINBistabilityModelCellNew")
        return new AuxinExoBistability(paraValue, indValue);
    else if (idValue == "SimpleROPModel")
        return new SimpleROPModel(paraValue, indValue);
    else if (idValue == "SimpleROPModel2")
        return new SimpleROPModel2(paraValue, indValue);
    else if (idValue == "SimpleROPModel3")
        return new SimpleROPModel3(paraValue, indValue);
    else if (idValue == "SimpleROPModel4")
        return new SimpleROPModel4(paraValue, indValue);
    else if (idValue == "SimpleROPModel5")
        return new SimpleROPModel5(paraValue, indValue);
    else if (idValue == "SimpleROPModel6")
        return new SimpleROPModel6(paraValue, indValue);
    else if (idValue == "SimpleROPModel7")
        return new SimpleROPModel7(paraValue, indValue);
    else if (idValue == "UpInternalGradientModel")
        return new UpInternalGradientModel(paraValue, indValue);
    else if (idValue == "DownInternalGradientModel")
        return new DownInternalGradientModel(paraValue, indValue);
    else if (idValue == "DownInternalGradientModelGeometric")
        return new DownInternalGradientModelGeometric(paraValue, indValue);
    else if (idValue == "DownInternalGradientModelSingleCell")
        return new DownInternalGradientModelSingleCell(paraValue, indValue);
    else if (idValue == "UpExternalGradientModel")
        return new UpExternalGradientModel(paraValue, indValue);
    else if (idValue == "UpInternalGradientModel")
        return new UpInternalGradientModel(paraValue, indValue);
    else if (idValue == "AuxinFluxModel")
        return new AuxinFluxModel(paraValue, indValue);
    else if (idValue == "IntracellularPartitioning")
        return new IntracellularPartitioning(paraValue, indValue);
    else if (idValue == "IntracellularCoupling")
        return new IntracellularCoupling(paraValue, indValue);
    else if (idValue == "IntracellularIndirectCoupling")
        return new IntracellularIndirectCoupling(paraValue, indValue);

    // directionReaction.h, directionUpdate.cc
    else if (idValue == "ContinousMTDirection")
        return new ContinousMTDirection(paraValue, indValue);
    else if (idValue == "ContinousMTDirection3d")
        return new ContinousMTDirection3d(paraValue, indValue);
    else if (idValue == "UpdateMTDirection")
        return new UpdateMTDirection(paraValue, indValue);
    else if (idValue == "UpdateMTDirectionEquilibrium")
        return new UpdateMTDirectionEquilibrium(paraValue, indValue);
    else if (idValue == "UpdateMTDirectionConcenHill")
        return new UpdateMTDirectionConcenHill(paraValue, indValue);
    else if (idValue == "RotatingDirection")
        return new RotatingDirection(paraValue, indValue);

    // sisterVertex.h, sisterVertex.cc
    else if (idValue == "SisterVertex::InitiateFromFile")
        return new SisterVertex::InitiateFromFile(paraValue, indValue);
    else if (idValue == "SisterVertex::InitiateFromDistance")
        return new SisterVertex::InitiateFromDistance(paraValue, indValue);
    else if (idValue == "SisterVertex::Spring")
        return new SisterVertex::Spring(paraValue, indValue);
    else if (idValue == "SisterVertex::SpringCellConc")
        return new SisterVertex::SpringCellConc(paraValue, indValue);
    else if (idValue == "SisterVertex::CombineDerivatives")
        return new SisterVertex::CombineDerivatives(paraValue, indValue);

    // calculate.h (.cc)
    // namespace Calculate collecting some ad hoc rections for calculating useful
    // information to be stored in cell or wall data
    else if (idValue == "Calculate::AngleVectors" ||
             idValue == "CalculateAngleVectors")
        return new Calculate::AngleVectors(paraValue, indValue);
    else if (idValue == "Calculate::AngleVectorXYplane" ||
             idValue == "CalculateAngleVectorXYplane")
        return new Calculate::AngleVectorXYplane(paraValue, indValue);
    else if (idValue == "Calculate::AngleVector" ||
             idValue == "AngleVector")
        return new Calculate::AngleVector(paraValue, indValue);
    else if (idValue == "Calculate::VertexVelocity")
        return new Calculate::VertexVelocity(paraValue, indValue);
    else if (idValue == "Calculate::TissueVolumeChange" ||
             idValue == "TemplateVolumeChange")
        return new Calculate::TissueVolumeChange(paraValue, indValue);

    // adhocReaction.h,adhocReaction.cc
    else if (idValue == "InflationDeflationStresses")
        return new InflationDeflationStresses(paraValue, indValue);
    else if (idValue == "VertexNoUpdateFromPosition")
        return new VertexNoUpdateFromPosition(paraValue, indValue);
    else if (idValue == "VertexNoUpdateFromIndex")
        return new VertexNoUpdateFromIndex(paraValue, indValue);
    else if (idValue == "VertexNoUpdateFromIndexHoldX")
        return new VertexNoUpdateFromIndexHoldX(paraValue, indValue);
    else if (idValue == "VertexNoUpdateFromIndexHoldY")
        return new VertexNoUpdateFromIndexHoldY(paraValue, indValue);
    else if (idValue == "VertexNoUpdateFromIndexHoldZ")
        return new VertexNoUpdateFromIndexHoldZ(paraValue, indValue);
    else if (idValue == "VertexNoUpdateFromList")
        return new VertexNoUpdateFromList(paraValue, indValue);
    else if (idValue == "VertexRandTip")
        return new VertexRandTip(paraValue, indValue);
    else if (idValue == "VertexNoUpdateBoundary")
        return new VertexNoUpdateBoundary(paraValue, indValue);
    else if (idValue == "VertexNoUpdateBoundaryPtemplate")
        return new VertexNoUpdateBoundaryPtemplate(paraValue, indValue);
    else if (idValue == "VertexNoUpdateBoundaryPtemplateStatic")
        return new VertexNoUpdateBoundaryPtemplateStatic(paraValue, indValue);
    else if (idValue == "VertexNoUpdateBoundaryPtemplateStatic3D")
        return new VertexNoUpdateBoundaryPtemplateStatic3D(paraValue, indValue);
    else if (idValue == "VertexNoUpdateBoundary3D")
        return new VertexNoUpdateBoundary3D(paraValue, indValue);
    else if (idValue == "VertexFromConstStressBoundary")
        return new VertexFromConstStressBoundary(paraValue, indValue);
    else if (idValue == "cellPolarity3D")
        return new cellPolarity3D(paraValue, indValue);
    else if (idValue == "diffusion3D")
        return new diffusion3D(paraValue, indValue);
    else if (idValue == "manipulate")
        return new manipulate(paraValue, indValue);
    else if (idValue == "VertexTranslateToMax")
        return new VertexTranslateToMax(paraValue, indValue);
    else if (idValue == "CenterCOM")
        return new CenterCOM(paraValue, indValue);
    else if (idValue == "CenterCOMcenterTriangulation")
      return new CenterCOMcenterTriangulation(paraValue, indValue);
    else if (idValue == "CenterCellCOM")
      return new CenterCellCOM(paraValue, indValue);
    else if (idValue == "CalculatePCAPlane")
        return new CalculatePCAPlane(paraValue, indValue);
    else if (idValue == "InitiateWallLength")
        return new InitiateWallLength(paraValue, indValue);
    else if (idValue == "InitiateTargetArea")
        return new InitiateTargetArea(paraValue, indValue);
    else if (idValue == "InitiateWallVariableConstant")
        return new InitiateWallVariableConstant(paraValue, indValue);
    else if (idValue == "InitiateWallMesh")
        return new InitiateWallMesh(paraValue, indValue);
    else if (idValue == "StrainTest")
        return new StrainTest(paraValue, indValue);
    else if (idValue == "CalculateVertexStressDirection")
        return new CalculateVertexStressDirection(paraValue, indValue);
    else if (idValue == "MoveVerticesRandomlyCapCylinder")
        return new MoveVerticesRandomlyCapCylinder(paraValue, indValue);
    else if (idValue == "scaleTemplate")
        return new scaleTemplate(paraValue, indValue);
    else if (idValue == "copyCellVector")
        return new copyCellVector(paraValue, indValue);
    else if (idValue == "randomizeMT")
        return new randomizeMT(paraValue, indValue);
    else if (idValue == "restrictVertexRadially")
        return new restrictVertexRadially(paraValue, indValue);
    else if (idValue == "CreationPrimordiaTime")
        return new CreationPrimordiaTime(paraValue, indValue);
    else if (idValue == "VertexFromRotationalForceLinear")
        return new VertexFromRotationalForceLinear(paraValue, indValue);
    else if (idValue == "ThresholdSwitch")
        return new ThresholdSwitch(paraValue, indValue);
    else if (idValue == "ThresholdReset")
        return new ThresholdReset(paraValue, indValue);
    else if (idValue == "ThresholdNoisyReset")
        return new ThresholdNoisyReset(paraValue, indValue);
    else if (idValue == "ThresholdResetAndCount")
        return new ThresholdResetAndCount(paraValue, indValue);
    else if (idValue == "FlagNoisyReset")
        return new FlagNoisyReset(paraValue, indValue);
    else if (idValue == "ThresholdAndFlagNoisyReset")
        return new ThresholdAndFlagNoisyReset(paraValue, indValue);
    else if (idValue == "FlagAddValue")
        return new FlagAddValue(paraValue, indValue);
    else if (idValue == "CopyVariable")
        return new CopyVariable(paraValue, indValue);
    else if (idValue == "DebugReaction")
        return new DebugReaction(paraValue, indValue);

    // Boolean namespace (boolean.h)
    else if (idValue == "Boolean::AndGate" ||
             idValue == "AndGate")
        return new Boolean::AndGate(paraValue, indValue);
    else if (idValue == "Boolean::AndNotGate" ||
             idValue == "AndNotGate")
        return new Boolean::AndNotGate(paraValue, indValue);
    else if (idValue == "Boolean::AndSpecialGate" ||
             idValue == "AndSpecialGate")
        return new Boolean::AndSpecialGate(paraValue, indValue);
    else if (idValue == "Boolean::AndSpecialGate2" ||
             idValue == "AndSpecialGate2")
        return new Boolean::AndSpecialGate2(paraValue, indValue);
    else if (idValue == "Boolean::AndSpecialGate3" ||
             idValue == "AndSpecialGate3")
        return new Boolean::AndSpecialGate3(paraValue, indValue);
    else if (idValue == "Boolean::AndThresholdGate" ||
             idValue == "AndThresholdsGate")
        return new Boolean::AndThresholdsGate(paraValue, indValue);
    else if (idValue == "Boolean::Count" ||
             idValue == "Count")
        return new Boolean::Count(paraValue, indValue);
    else if (idValue == "Boolean::FlagCount" ||
             idValue == "FlagCount")
        return new Boolean::FlagCount(paraValue, indValue);
    else if (idValue == "Boolean::AndGateCount" ||
             idValue == "AndGateCount")
        return new Boolean::AndGateCount(paraValue, indValue);
    else if (idValue == "Boolean::OrGateCount" ||
             idValue == "OrGateCount")
        return new Boolean::OrGateCount(paraValue, indValue);
    else if (idValue == "Boolean::OrSpecialGateCount" ||
             idValue == "OrSpecialGateCount")
        return new Boolean::OrSpecialGateCount(paraValue, indValue);

    // cellTime.h
    else if (idValue == "CellTimeDerivative")
        return new CellTimeDerivative(paraValue, indValue);

    // MembraneCycling.h
    else if (idValue == "MembraneCycling::Constant")
        return new MembraneCycling::Constant(paraValue, indValue);
    else if (idValue == "MembraneCycling::CrossMembraneNonLinear")
        return new MembraneCycling::CrossMembraneNonLinear(paraValue, indValue);
    else if (idValue == "MembraneCycling::LocalWallFeedbackNonLinear")
        return new MembraneCycling::LocalWallFeedbackNonLinear(paraValue, indValue);
    else if (idValue == "MembraneCycling::LocalWallFeedbackLinear")
        return new MembraneCycling::LocalWallFeedbackLinear(paraValue, indValue);
    else if (idValue == "MembraneCycling::CellUpTheGradientNonLinear")
        return new MembraneCycling::CellUpTheGradientNonLinear(paraValue, indValue);
    else if (idValue == "MembraneCycling::CellUpTheGradientLinear")
        return new MembraneCycling::CellUpTheGradientLinear(paraValue, indValue);
    else if (idValue == "MembraneCycling::PINFeedbackNonLinear")
        return new MembraneCycling::PINFeedbackNonLinear(paraValue, indValue);
    else if (idValue == "MembraneCycling::PINFeedbackLinear")
        return new MembraneCycling::PINFeedbackLinear(paraValue, indValue);
    else if (idValue == "MembraneCycling::InternalCellNonLinear")
        return new MembraneCycling::InternalCellNonLinear(paraValue, indValue);
    else if (idValue == "MembraneCycling::InternalCellLinear")
        return new MembraneCycling::InternalCellLinear(paraValue, indValue);
    else if (idValue == "MembraneCycling::CellFluxExocytosis")
        return new MembraneCycling::CellFluxExocytosis(paraValue, indValue);

    // MembraneCyclingAll.h
    else if (idValue == "MembraneCyclingAll::Constant")
        return new MembraneCyclingAll::Constant(paraValue, indValue);
    else if (idValue == "MembraneCyclingAll::LocalWallFeedbackNonLinear")
        return new MembraneCyclingAll::LocalWallFeedbackNonLinear(paraValue, indValue);
    else if (idValue == "MembraneCyclingAll::LocalWallFeedbackNonLinearInhibition")
        return new MembraneCyclingAll::LocalWallFeedbackNonLinearInhibition(paraValue, indValue);

    // massAction.h
    else if (idValue == "MassAction::General")
        return new MassAction::General(paraValue, indValue);
    else if (idValue == "MassAction::GeneralEnzymatic")
        return new MassAction::GeneralEnzymatic(paraValue, indValue);
    else if (idValue == "MassAction::OneToTwo")
        return new MassAction::OneToTwo(paraValue, indValue);
    else if (idValue == "MassAction::TwoToOne")
        return new MassAction::TwoToOne(paraValue, indValue);
    else if (idValue == "MassAction::GeneralWall")
        return new MassAction::GeneralWall(paraValue, indValue);
    else if (idValue == "MassAction::OneToTwoWall")
        return new MassAction::OneToTwoWall(paraValue, indValue);
    else if (idValue == "MassAction::TwoToOneWall")
        return new MassAction::TwoToOneWall(paraValue, indValue);
    else if (idValue == "MassAction::HillSimple")
        return new MassAction::HillSimple(paraValue, indValue);

    // Namespace Initiation, initiation.h
    else if (idValue == "Initiation::Random")
        return new Initiation::Random(paraValue, indValue);
    else if (idValue == "Initiation::RandomBoolean")
        return new Initiation::RandomBoolean(paraValue, indValue);
    else if (idValue == "Initiation::RandomBooleanBiased")
        return new Initiation::RandomBooleanBiased(paraValue, indValue);
    else if (idValue == "Initiation::FaceArea2D")
        return new Initiation::FaceArea2D(paraValue, indValue);
    else if (idValue == "Initiation::CellLabel")
        return new Initiation::CellLabel(paraValue, indValue);

    // Namespace Hypocotyl3D, hypocotyl3D.h
    else if (idValue == "Hypocotyl3D::limitZdis")
        return new Hypocotyl3D::limitZdis(paraValue, indValue);
    else if (idValue == "Hypocotyl3D::StrainTRBS")
        return new Hypocotyl3D::StrainTRBS(paraValue, indValue);
    else if (idValue == "Hypocotyl3D::VertexFromTRBScenterTriangulationMT")
        return new Hypocotyl3D::VertexFromTRBScenterTriangulationMT(paraValue, indValue);

    // Obselete reactions
    if (idValue == "MembraneDiffusionSimple2") {
        std::cerr << "BaseReaction::createReaction() EXITING: "
                  << "Reaction MembraneDiffusionSimple2 has been removed (20210820) since only working on"
                  << " a single cell."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (idValue == "Calculate::MaxVelocity" ||
        idValue == "maxVelocity") {
        std::cerr << "BaseReaction::createReaction() EXITING: "
                  << "Reaction Calculate::MaxVelocity (maxVelocity) has been replaced by Calculate::VertexVelocity."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (idValue == "WallGrowthExponentialTruncated") {
        std::cerr << "BaseReaction::createReaction() EXITING: "
                  << "Reaction WallGrowthExponentialTruncated has been replaced by WallGrowth::Constant."
                  << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "VertexFromHypocotylGrowth") {
        std::cerr << "BaseReaction::createReaction() EXITING: "
                  << "Reaction VertexFromHypocotylGrowth "
                  << "has been replaced by Force::Axial. "
                  << "See documentation for alternative versions."
                  << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "EpidermalVertexForce") {
        std::cerr << "BaseReaction::createReaction() EXITING: "
                  << "Reaction EpidermalVertexForce "
                  << "has been replaced by Force::EpidermalCoordinate."
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    else if (idValue == "WallGrowthExponentialStressTruncated") {
        std::cerr << "BaseReaction::createReaction() EXITING: "
                  << "Reaction WallGrowthExponentialStressTruncated "
                  << "has been replaced by WallGrowth::Stress (setting "
                  << "the stretch_flag to 1 and provide L_th)."
                  << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "WallLengthGrowExperimental") {
        std::cerr << "BaseReaction::createReaction() EXITING: "
                  << "Reaction WallLengthGrowExperimental "
                  << "has been replaced by WallGrowth::Force. Better is to use the WallGrowth::Stress "
                  << "(setting the stretch_flag to 0 and not provide L_th)."
                  << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "WallGrowthConstantStress" ||
               idValue == "WallGrowthConstantStressConcentrationHill") {
        std::cerr << "BaseReaction::createReaction() EXITING: "
                  << "WallGrowthConstantStress* has been "
                  << "replaced by WallGrowth::Stress*." << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "VertexFromWallSpringExperimental") {
        std::cerr << std::endl
                  << "BaseReaction::createReaction() EXITING: "
                  << "Reaction VertexFromWallSpringExperimental not used anymore." << std::endl
                  << "Use WallMechanics::Spring with parameter(2) set to 1.0 instead." << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "VertexFromWallSpringAsymmetric" ||
               idValue == "VertexFromEpidermalWallSpringAsymmetric" ||
               idValue == "VertexFromEpidermalCellWallSpringAsymmetric") {
        std::cerr << "BaseReaction::BaseReaction() EXITING: "
                  << "All reactions named *SpringAsymmetric have been renamed "
                  << "*Spring." << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "VertexFromCellPressure") {
        std::cerr << std::endl
                  << "BaseReaction::createReaction() EXITING: "
                  << "Reaction VertexFromCellPressure renamed." << std::endl
                  << "Use Pressure2D::EdgeForce." << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "VertexFromCellPressureVolumeNormalized") {
        std::cerr << std::endl
                  << "BaseReaction::createReaction() EXITING: "
                  << "Reaction VertexFromCellPressureVolumeNormalized renamed (and reimplemented)." << std::endl
                  << "Use Pressure2D::AreaPotential with parameter(1) set to 1 instead." << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "VertexFromCellPressureThresholdFromMaxPos") {
        std::cerr << std::endl
                  << "BaseReaction::createReaction() EXITING: "
                  << "Reaction VertexFromCellPressureThresholdFromMaxPos renamed." << std::endl
                  << "Use Pressure2D::AreaPotentialSpatialThreshold instead." << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "VertexFromCellInternalPressure") {
        std::cerr << std::endl
                  << "BaseReaction::createReaction() EXITING: "
                  << "Reaction VertexFromCellInternalPressure not used anymore." << std::endl
                  << "Use Pressure2D::AreaPotential with parameter(2) set to 1 instead." << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "VertexFromPressureExperimental") {
        std::cerr << std::endl
                  << "BaseReaction::createReaction() EXITING: "
                  << "Reaction VertexFromPressureExperimental renamed to "
                  << "Pressure2D::AreaPotentialTargetArea." << std::endl;
        exit(EXIT_FAILURE);
    } else if (idValue == "CellVolumeExperimental") {
        std::cerr << std::endl
                  << "BaseReaction::createReaction() EXITING: "
                  << "Reaction CellVolumeExperimental renamed to "
                  << "TargetAreaFromPressure." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Default, if nothing found
    else {
        std::cerr << std::endl
                  << "BaseReaction::createReaction() EXITING: "
                  << "Reactiontype "
                  << idValue << " not known, no reaction created." << std::endl;
        exit(EXIT_FAILURE);
    }
}

BaseReaction *
BaseReaction::createReaction(std::istream &IN) {
    std::string idVal;
    size_t pNum, levelNum;
    IN >> idVal;
    IN >> pNum;
    IN >> levelNum;
    std::vector<size_t> varIndexNum(levelNum);
    for (size_t i = 0; i < levelNum; i++)
        IN >> varIndexNum[i];

    std::vector<double> pVal(pNum);
    for (size_t i = 0; i < pNum; i++)
        IN >> pVal[i];

    std::vector<std::vector<size_t> > varIndexVal(levelNum);
    for (size_t i = 0; i < levelNum; i++)
        varIndexVal[i].resize(varIndexNum[i]);

    for (size_t i = 0; i < levelNum; i++)
        for (size_t j = 0; j < varIndexNum[i]; j++)
            IN >> varIndexVal[i][j];

    return createReaction(pVal, varIndexVal, idVal);
}

void BaseReaction::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &walldata,
           DataMatrix &vertexData,
           DataMatrix &cellderivs,
           DataMatrix &wallderivs,
           DataMatrix &vertexDerivs) {
    std::cerr << "BaseReaction::derivs() should not be used. "
              << "Should always be mapped onto one of the real types." << std::endl;
    exit(0);
}

void BaseReaction::
    derivsWithAbs(Tissue &T,
                  DataMatrix &cellData,
                  DataMatrix &walldata,
                  DataMatrix &vertexData,
                  DataMatrix &cellderivs,
                  DataMatrix &wallderivs,
                  DataMatrix &vertexDerivs,
                  DataMatrix &sdydtCell,
                  DataMatrix &sdydtWall,
                  DataMatrix &sdydtVertex) {
    std::cerr << "BaseReaction::derivsWithAbs() should not be used. "
              << "Should always be mapped onto one of the real types." << std::endl;
    exit(0);
}

void BaseReaction::initiate(Tissue &T,
                            DataMatrix &cellData,
                            DataMatrix &walldata,
                            DataMatrix &vertexData,
                            DataMatrix &cellderivs,
                            DataMatrix &wallderivs,
                            DataMatrix &vertexDerivs) {
}

void BaseReaction::update(Tissue &T,
                          DataMatrix &cellData,
                          DataMatrix &walldata,
                          DataMatrix &vertexData,
                          double h) {
}

void BaseReaction::print(std::ofstream &os) {
    std::cerr << "BaseReaction::print(ofstream) should not be used. "
              << "Should always be mapped onto one of the real types.\n";
    exit(EXIT_FAILURE);
}

void BaseReaction::printPly(std::ofstream &os) {
    // Only implemented for some reactions, do nothing for others
}

void BaseReaction::printState(Tissue *T,
                              DataMatrix &cellData,
                              DataMatrix &wallData,
                              DataMatrix &vertexData,
                              std::ostream &os) {
}
