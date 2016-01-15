// System Includes
#include <iostream>
#include <sstream>
#include <string>
#include "ReportRipperStandalone.hpp"


#ifdef  __CLASS__
#undef  __CLASS__
#endif
#define __CLASS__ "main"


// Namespaces used
using namespace std;


int main(int argc, char ** argv)
{
  string file;
  if (argc == 2) 
  {
    file = string(argv[1]);
    if (file.at(0) == '-' && file.size() > 1) 
    {
      cerr << "Usage: stitch [m8 file]" << endl;
      exit(-1);
    }
  }

  ReportRipperStandalone rrs(file);

  cout << rrs.print() << flush;
}
