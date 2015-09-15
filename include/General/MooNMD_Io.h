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
#include <stdexcept>

using std::cout;
using std::cerr;
using std::endl;
using std::ends;
using std::setw;
using std::setprecision;

void OpenFiles();
void CloseFiles();

#define OutPut(x) {OutFile << x; cout << x;}
#define ErrMsg(x) { cerr <<"Error in file "<<__FILE__<<", line "<<__LINE__<<":\n"<<x<<endl; OutFile<<"Error in file "<<__FILE__<<", line "<<__LINE__<<":\n"<<x<<endl;}
#define Error(x) {ErrMsg(x)}
void throw_with_message(const std::string& x, std::string file, int line);
#define ErrThrow(x) {throw_with_message(std::string() + x, __FILE__, __LINE__);}

extern std::ofstream OutFile;

#endif
