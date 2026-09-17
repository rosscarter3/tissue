//
// Filename     : force.h
// Description  : Classes describing mechanical updates from external forces
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : April 2018
// Revision     : $Id:$
//
#ifndef FORCE_H
#define FORCE_H

#include <cmath>

#include "baseReaction.h"
#include "tissue.h"

///
/// @brief Reactions describing updates coming from forces applied in different scenarios.
///
/// @details External forces can be applied to the tissue to represent walls (not allowed to pass)
/// or attraction to specific geometries to help simulations generate specific shapes etc.
/// These reaction are collected within this namespace.
///
namespace Force {
///
/// @brief Applies a force on each vertex proportianal to the distance of that vertex from a given
/// plane and perpendicular to that plane. For example, if y force required a force will be
/// generated in the y direction proportional to distance from xz plane in dircection
///
/// @details THis reaction applies a force to each vertex in the tissue proportional to the vertices
/// distance from a given axis (only really makes sense for a template where the centroid is at
/// (0,0,0)). This force is in the given direction perpendicualr to the plane formed by the other two
/// directions.
///
/// @f[ \frac{dx(y,z)}{dt} = p_0 * x(y,z) @f]
///
/// In a model file it is defined as: (UPDATE)
/// @verbatim
/// Force::ForceFromPlane 2 0
/// F
/// direction_flag (0:x, 1:y, 2:z)
/// @endverbatim
/// where F is the force magnitude and the direction flag is the required force direction (1,2,3)/(x,y,z)
///
///
///
class ForceFromPlane : public BaseReaction {
   public:
    ForceFromPlane(std::vector<double> &paraValue,
                   std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Applies a force outwards (or inwards) perpendicular to a Cylinder surface direction
///
/// @details This reaction applies a force in the radial (x,y) direction with
/// update of the vertex position in x (and symmetric in y) following
/// @f[ \frac{dx}{dt} = p_0 p_1 \frac{x}{|r|} @f]
/// where @f$p_0@f$ is the force and @f$p_1@f$ sets the direction (1 outwards, -1 inwards).
/// @f$|r|@f$ is the radial position (in the x/y-plane).
/// In a model file it is defined as:
/// @verbatim
/// Force::Cylinder 2 0
/// F direction_flag
/// @endverbatim
/// where F is the force magnitude and the direction_flag is 1 if forces are radially
/// (x,y) outwards and -1 for inwards.
///
/// @note Used to be called CylinderForce
///
class Cylinder : public BaseReaction {
   public:
    Cylinder(std::vector<double> &paraValue,
             std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Applies a force outwards (or inwards) perpendicular to a SphereCylinder surface direction
///
/// @details This reaction applies a force in the radial direction, where z>0 updates vertex positions
/// in x,y,z (sphere), and z<0 only x,y directions are updated (cylinder).
/// For z>0, the update of the vertex position in x (and symmetric in y and z) follows
/// @f[ \frac{dx}{dt} = p_0 p_1 \frac{x}{|r|} @f]
/// where @f$p_0@f$ is the force and @f$p_1@f$ sets the direction (1 outwards, -1 inwards).
/// @f$|r|@f$ is the radial position (in 3D (x,y,z)).
/// For z<0, the update of the vertex position in x (and symmetric in y) follows
/// @f[ \frac{dx}{dt} = p_0 p_1 \frac{x}{|r|} @f]
/// where @f$p_0@f$ is the force and @f$p_1@f$ sets the direction (1 outwards, -1 inwards).
/// @f$|r|@f$ is now the radial position in the x/y-plane.
/// In a model file it is defined as:
/// @verbatim
/// Force::SphereCylinder 2 0
/// F direction_flag
/// @endverbatim
/// where F is the force magnitude and the direction_flag is 1 if forces are radially
/// (x,y) outwards and -1 for inwards.
///
/// @note Used to be called SphereCylinderForce
/// @see GrowthForce::SphereCylinder for generating vertex movement alonf the sphere-cylinder surface.
///
class SphereCylinder : public BaseReaction {
   public:
    SphereCylinder(std::vector<double> &paraValue,
                   std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Applies a force perpendicular to a spherecylinder surface, outwards (inwards) if vertex is
/// inside (outside) a specified radius.
///
/// @details This reaction generates forces perpendicular to a sphere (if z>0) or cylinder (z<0) surface
/// of specified radius (@f$R=p_2@f$).
/// For z>0, the update of the vertex position in x (and symmetric in y and z) follows
/// (if vertex inside of the radius)
/// @f[ \frac{dx}{dt} = p_0 \frac{p_2-|r|}{|r|} x @f]
/// where @f$p_0@f$ is the force and @f$p_2@f$ is the defined radius. If outside @f$p_1@f$ replaces @f$p_0@f$
/// and direction is inwards. @f$|r|@f$ is the radial position of the vertex (in 3D (x,y,z)).
/// For z<0, the update of the vertex position in x (and symmetric in y) follows
/// @f[ \frac{dx}{dt} = p_0 \frac{p_2-|r|}{|r|} x @f]
/// with same parameters and on the sphere, but
/// @f$|r|@f$ is now the radial position in the x/y-plane.
/// In a model file it is defined as:
/// @verbatim
/// Force::SphereCylinderRadius 3 0
/// F_out F_in R
/// @endverbatim
/// where F_out is the force magnitude applied (outwards) if the vertex is inside the specified radius (R),
/// and F_in is if it is outside (in inwards direction).
///
/// @note Used to be called SphereCylinderForceFromRadius
/// @see GrowthForce::SphereCylinder applies forces along a spherecylinder surface.
///
class SphereCylinderRadius : public BaseReaction {
   public:
    SphereCylinderRadius(std::vector<double> &paraValue,
                         std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Applies a force perpendicular to a defined wall of infinite size, defined for a specific coordinate
///
/// @details A spring force in a perpendicular direction to a specific coordinate is applied.
/// The update is given by
/// @f[\frac{dx[/y/z]_i}{dt} -= p_0 (x_i-p_1) @f]
/// if @f$x_i>p_1@f$, i.e. the vertex has crossed the 'wall', and for the coordinate specified.
/// This is when @f$p_2=1@f$. When @f$p_2=-1@f$, the update is
/// @f[\frac{dx[/y/z]_i}{dt} -= p_0 (x_i-p_1) @f]
/// if @f$x_i<p_1@f$, i.e. crossing the 'wall' from the right.
/// In a model file, the reaction is given by
/// @verbatim
/// Force::InfiniteWall 3 1 1
/// K_force threshold(wall position) direction_flag
/// position(coordinate)
/// @endverbatim
///
/// @see Force::ExternalWall for more complex version of a wall force that can be defined
/// in more general direction and that can move.
class InfiniteWall : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    InfiniteWall(std::vector<double> &paraValue,
                 std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Applies a force on epidermal vertices in direction (coordinate) given by index
///
/// @details A force in a specific coordinate direction is applied to epidermal vertices.
/// It will update epidermal indices according to
/// @f[\frac{dx[y,z]_i}{dt} += p_0 p_1 @f]
/// where the coordinate is given as first index and @f$p_0@f$ is the force and @f$p_1@f$
/// is a flag +/-1 to set the direction (1->outwards, -1->inwards).
/// In a model file this is given by
/// @verbatim
/// Force::EpidermalCoordinate 2 1 1
/// K_force direction(=1/-1)
/// coordinate_index
/// @endverbatim
///
class EpidermalCoordinate : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    EpidermalCoordinate(std::vector<double> &paraValue,
                        std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Applies a force on epidermal vertices in a radial direction
///
/// @details A force in a radial direction is applied to epidermal vertices.
/// It will update epidermal indices according to
/// @f[\frac{dx[y]_i}{dt} += p_0 \frac{x_i}{R_i} @f]
/// where x and y coordinates for vertex i are updated. @f$p_0@f$ is the force and @f$R_i@f$
/// is the radius (distance to origo (x,y)) of the vertex.
/// It moves vertices outwards if the parameter is positive (and inwards if negative).
/// In a model file this is given by
/// @verbatim
/// Force::EpidermalRadial 1 0
/// K_force
/// @endverbatim
///
/// @see Force::Radial
///
class EpidermalRadial : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    EpidermalRadial(std::vector<double> &paraValue,
                    std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Applies a force in a radial direction on vertices listed
///
/// @details A force in a radial direction is applied to vertices listed.
/// It will update epidermal indices according to
/// @f[\frac{dx[y,z]_i}{dt} += p_0 p_1 x_i @f]
/// where each coordinate is updated @f$p_0@f$ is the force.
/// Note, @f$p_1@f$ is a flag setting the direction i.e
/// if 1, the reaction moves vertices outwards (and inwards if -1).
/// In a model file this is given by
/// @verbatim
/// Force::IndexRadial 1 1 N
/// K_force
/// v_1 ... v_N
/// @endverbatim
///
/// @note Used to be named VertexForceOrigoFromIndex.
/// @see Force::Radial
/// @see Force::EpidermalRadial (similar update but vertices are selected to be updated if they are at boundary)
/// @see Force::CellIndexRadial (similar but cell indices are listed to set vertices to be updated)
///
class IndexRadial : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    IndexRadial(std::vector<double> &paraValue,
                std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Applies a force in a radial direction on vertices of cells listed
///
/// @details A force in a radial direction is applied to vertices of cells listed.
/// It will update epidermal indices according to
/// @f[\frac{dx[y,z]_i}{dt} += p_0 p_1 x_i @f]
/// where each coordinate is updated @f$p_0@f$ is the force.
/// Note, @f$p_1@f$ is a flag setting the direction i.e
/// if 1, the reaction moves vertices outwards (and inwards if -1).
/// In a model file this is given by
/// @verbatim
/// Force::CellIndexRadial 1 1 N
/// K_force direction_flag(=+/-1)
/// c_1 ... c_N
/// @endverbatim
/// where the list of cells indices whose vertices should be updated.
///
/// @note Used to be named CellForceOrigoFromIndex.
/// @see Force::Radial
/// @see Force::EpidermalRadial (similar update but vertices are selected to be updated if they are at boundary)
/// @see Force::CellIndexRadial (similar but cell indices are listed to set vertices to be updated)
///
class CellIndexRadial : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    CellIndexRadial(std::vector<double> &paraValue,
                    std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Applies a constant axial force on vertices in two intervals in z (opposing force direction),
/// used e.g on a cylinderical template to approximate internal tissue forces generating axial forces/growth.
///
/// @details A constant force, F, is applied axially to regions in z in [z_0+a, z_0+a+d] (upward force) and
/// [z_0-a, z_0-a-d] (downward) resembling forces generated by internal tissues in e.g.
/// Hypocotyls that can provide forces for axial growth.
/// In a model file the reaction is defined as:
/// @verbatim
/// Force::Axial 4 0
/// z_0 a d F
/// @endverbatim
/// where z_0 is the 'center' of the tissue, a and d defines the two regions, and F is the magnitude of the F.
///
/// @note Used to be called VertexFromHypocotylGrowth
/// @note For developers: A more advanced version (8 parmeters) generating axial and radial forces is commented out in force.h(.cc)
///
class Axial : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    Axial(std::vector<double> &paraValue,
          std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Updates list of vertices with a given force in direction provided by vector
///
/// @details A force is given as parameters (p[0] for Fx, [optional p[1/2] for Fy/Fz].
/// The update is given by (for dimensions i as specified by the parameters)
/// @f[ \frac{dx_i}{dt} = p_i @f]
/// In a model file the reaction is defined as
/// @verbatim
/// VertexFromForce 1[/2/3] 1 num_vertices
/// F_x [F_y F_z]
/// vertex_index_0
/// [vertex index_1]
/// [...]
/// @endverbatim
/// @note Used to be called VertexFromForce
/// @see Force::VectorLinear Same as this reaction, but uses a 'ramping up' of the forces (linearly) over dT.
///
class Vector : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    Vector(std::vector<double> &paraValue,
           std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Updates list of vertices with a given force applied where the force
/// is linearly increased from zero across a given time span (deltaT).
///
/// The update is given by (for dimensions i as specified by the parameters)
/// @f[ \frac{dx_i}{dt} = max{\frac{t}{T},1.0} p_i @f]
/// where @f$p_i@f$ are the forces and the factor in front is between 0 and 1
/// linear increasing between time (t) equals zero and t equals T.
/// In a model file the reaction is defined as
/// @verbatim
/// Force::VectorLinear 2[3/4] 1 num_vertices
/// F_x [F_y F_z] T
/// vertex_index_0
/// [vertex index_1]
/// [...]
/// @endverbatim
///
/// @note Used to be called VertexFromForceLinear.
/// @see Force::Vector Same as this reaction, but with Forces constant.
///
class VectorLinear : public BaseReaction {
   private:
    double timeFactor_;

   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    VectorLinear(std::vector<double> &paraValue,
                 std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

///
/// @brief Updates position of vertices assuming that a ball (force) pushing at the tissue,
/// where it is possible to have the ball moving at a given constant velocity
///
/// @details The force is applied outwards with respect to a sphere/ball and the force
/// is proportional to
/// (overlap)^(3/2). In a model file the reaction is defined as
/// @verbatim
/// Force::Ball 5 0
/// Radius Xc Yc Zc Kforce
/// @endverbatim
/// or
/// @verbatim
/// Force::Ball 8 0
/// Radius Xc Yc Zc Kforce dXc dYc dZc
/// @endverbatim
/// where radius is the size of the 'ball' pushing at the tissue, Xc,Yc,Zc is
/// the center of the ball, and the optional dXc,dYc,dZc are the rates for
/// moving the ball along the different directions (the movement is generated in
/// the update function).
///
/// @note Used to be called VertexFromBall
///
class Ball : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    Ball(std::vector<double> &paraValue,
         std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

///
/// @brief Updates position of vertices assuming that a parabolid object is pushed into
/// a tissue, and can be moving with a given velocity (in z) into the template.
///
/// @details The force is applied outwards with respect to the parabolid.
/// In a model file the reaction is defined as
/// @verbatim
/// Force::Parabolid 5 0
/// a Xc Yc b Kforce
/// @endverbatim
/// or
/// @verbatim
/// Force::Parabolid 8 0
/// a Xc Yc b Kforce VelocityZ
/// @endverbatim
/// where the parabolid is defined by z=a((x-Xc)2 +(y-Yc)2)+b
/// Xc,Yc is the center in x,y of the parabolid, and a and b are size/shape parameters.
/// The optional VelocityZ is the rate for moving the ball along the z axis.
/// (the movement is generated in the update function).
///
/// @note Used to be called VertexFromParabolid.
///
class Parabolid : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    Parabolid(std::vector<double> &paraValue,
              std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

///
/// @brief Updates position of vertices assuming that an external wall is moving
/// with a given velocity vector toward the meristem.
///
/// @details The force is applied in the normal direction and is proportional to
/// (overlap)^(3/2) where the force magnitude is given by p10. (p0, p1, p2) defines the
/// original point on the plane, (p3, p4, p5) defines the normal of the plane, (p6, p7, p8)
/// defines the movement vector, and p9 defines the maximal distance of movement (measured
/// as the distance between current plane position point from original, which should be the
/// same as (sum dt)*|dX|).
/// In a model file the reaction is defined as
/// @verbatim
/// Force::ExternalWall 12 0
/// X0 Y0 Z0 nx ny nz dXc dYc dZc Dmax Kforce
/// @endverbatim
/// where n is the normal vector to the 'wall' pushing at the tissue, X0,Y0,Z0
/// is a point on the wall and dXc,dYc,dZc are the rates for moving the wall
/// along the different directions (the movement is generated in the update
/// function).
///
/// @note used to be called VertexFromExternalWall, and have a Zmin/max instead of Dmax
/// @see Force::InfiniteWall is a simpler version with fewer parameters, but with wall along
/// axis and static
///
class ExternalWall : public BaseReaction {
   private:
    // variables to store original positions (for calculation of plane distance movement)
    double x00_;
    double y00_;
    double z00_;

   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    ExternalWall(std::vector<double> &paraValue,
                 std::vector<std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
    ///
    /// @brief Prints the wall plane into a vtu file
    ///
    void printPly(std::ofstream &os);
};

}  // end namespace Force

#endif
