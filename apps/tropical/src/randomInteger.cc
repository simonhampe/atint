#include "polymake/client.h"
#include "polymake/Array.h"
#include "polymake/RandomGenerators.h"

namespace polymake { namespace tropical {

  Array<Integer> randomInteger(const int& max_arg, const int &n) {
    static Integer upperBound = 0;
    static UniformlyRandomRanged<Integer> rg(max_arg);
      if(max_arg != upperBound)  {
	rg = UniformlyRandomRanged<Integer>(max_arg);
	upperBound = max_arg;
      }
    Array<Integer> result(n);
    for(int i = 0; i < n; i++) {
	result[i] = rg.get();
    }
    return result;
  }
  
  
  
  UserFunction4perl("# @category Random number generators"
                  "# Returns n random integers in the range 0.. (max_arg-1),inclusive"
                  "# Note that this algorithm is not optimal for real randomness:"
                  "# If you change the range parameter and then change it back, you will"
                  "# usually get the exact same sequence"
                  "# @param Integer max_arg The upper bound for the random integers"
                  "# @param Integer n The number of integers to be created"
                  "# @return Array<Integer>",
                  &randomInteger,"randomInteger(Integer, Integer)");                  
  
}}