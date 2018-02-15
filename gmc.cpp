#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <utility>
#include <iostream>

#include "util/global.h"
#include "simulation/Simulation.h"

using namespace std;
using namespace Util;
using namespace GridMC;

/**
* Main program for GridMC summation.
*/
int main(int argc, char **argv)
{

   if (argc < 2) {
      cout << "Parameter file not provided; do nothing." << endl;
   } else {
      Simulation simulation(argv[1]);
      simulation.run();
   }

   return 0;
}
