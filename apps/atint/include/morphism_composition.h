/*
 T his program is free s*oftware; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor,
 Boston, MA  02110-1301, USA.
 
 ---
 Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 */

#ifndef MORPHISM_COMPOSITION_H
#define MORPHISM_COMPOSITION_H

namespace polymake { namespace atint { 
  
  /**
    @brief Computes the composition g(f) of two morphisms f and g (as in f:X->Y, g:Y->Z). Actually, f and g can also be piecewise linear maps on their domains, the method will work equally well. It naturally requires that the image of f is contained in the domain of g (except if f is global and surjective)
    @param perl::Object f A Morphism
    @param perl::Object g A Morphism, whose DOMAIN contains the image of f - except if f is a global affine linear and surjective map. In this case, g's DOMAIN needn't be the whole space. The composition will then be defined on the preimage of g's DOMAIN
    @return perl::Object A Morphism object, the composition "g after f", always defined on a homogenous domain
   */
   perl::Object morphism_composition(perl::Object f, perl::Object g);
  
}}

#endif // MORPHISM_COMPOSITION_H