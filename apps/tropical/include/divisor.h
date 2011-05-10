#include "polymake/client.h"

#ifndef ATINT_DIVISOR_H
#define ATINT_DIVISOR_H

/**
@brief Takes two fans and computes the intersection of both. The function relies on the fact that the latter fan is complete (i.e. its support is the whole ambient space) to compute the intersection correctly.
@param fan An arbitrary polyhedral fan
@param completeFan A complete polyhedral fan
@return fan::PolyhedralFan The intersection of both fans (whose support is equal to the support of fan). The 
resulting fan uses homogeneous coordinates if and only fan does. If fan has a property TROPICAL_WEIGHTS, 
the tropical weights of the refinement are also computed. If fan is zero-dimensional (i.e. a point), fan is returned.
*/
perl::Object intersect_complete_fan(perl::Object fan, perl::Object completeFan);

#endif