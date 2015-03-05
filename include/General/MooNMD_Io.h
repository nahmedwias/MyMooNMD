// =======================================================================
// 
// Purpose:     Collection of routine and methods for IO
//
// Author:      Gunar Matthies
//
// Date:        2004/02/19
// 
// =======================================================================

#ifndef __MOONMD_IO__
#define __MOONMD_IO__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using std::cout;
using std::cerr;
using std::endl;
using std::ends;
using std::setw;
using std::setprecision;

void OpenFiles();
void CloseFiles();

#define OutPut(x) {OutFile << x; cout << x;}
#define Error(x) {LogFile << x; cerr << x;}

extern std::ofstream OutFile;
extern std::ofstream LogFile;

#endif
