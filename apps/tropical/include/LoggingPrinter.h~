#include "polymake/client.h"

#ifndef ATINT_LOGGING_PRINTER_H
#define ATINT_LOGGING_PRINTER_H

/**
 @brief This class defines a logger that can be turned on or off in the following way: Put "using namespace xxx;" in your code and replace xxx with atint::donotlog, atint::dolog, atint::dotrace, depending on what you want. Under donotlog, all  calls to dbglog or dbgtrace have no effect. Unter dolog, only dbglog produces output and under dbgtrace, both commands produce output (to pm::cout). The commands are used as cout, i.e. dbglog << ... and dbgtrace << ...
*/

namespace polymake { namespace tropical {

/**
 @brief A dummy streambuffer class that does absolutely nothing. Used for deactivating logging.
*/
class DummyBuffer : public std::streambuf{
    
};

static DummyBuffer dummybf;
static std::ostream dbgstream(&dummybf);


}}

namespace atint {
  
  namespace donotlog{
    /**
      @brief A logger printer that outputs nothing
    */
    static pm::PlainPrinter<> dbglog(polymake::tropical::dbgstream); 
    /**
      @brief A logger printer that outputs nothing
    */
    static pm::PlainPrinter<> dbgtrace(polymake::tropical::dbgstream);
  }
  namespace dolog{
    /**
      @brief A logger printer, that prints to pm::cout
    */
    static pm::PlainPrinter<> dbglog(pm::cout);
    /**
      @brief A logger printer that outputs nothing
    */
    static pm::PlainPrinter<> dbgtrace(polymake::tropical::dbgstream);
  }
  namespace dotrace{
    /**
      @brief A logger printer, that prints to pm::cout
    */
    static pm::PlainPrinter<> dbglog(pm::cout);
    /**
      @brief A logger printer that prints to pm::cout
    */
    static pm::PlainPrinter<> dbgtrace(pm::cout);
  }
}
#endif