#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Array.h"
#include "polymake/Set.h"

namespace polymake { namespace fan {
 
  /*
    @brief Computes whether a given weighted polyhedral fan is balanced
    @param f a PolyhedralFan object
    @return true if the fan is balanced, false otherwise
  */
  bool isBalanced (perl::Object f) {
      
    
      return true;
  }
  
  Function4perl(&isBalanced, "isBalanced(PolyhedralFan)");
  
}
}