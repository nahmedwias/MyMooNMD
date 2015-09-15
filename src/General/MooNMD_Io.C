#include <Database.h>
#include <MooNMD_Io.h>

std::ofstream OutFile;

void OpenFiles()
{
  OutFile.open(TDatabase::ParamDB->OUTFILE);
  OutFile.setf(std::ios::scientific);
}

void CloseFiles()
{
  OutFile.close();
}

void throw_with_message(const std::string& x, std::string file, int line)
{
  std::string err_string = std::string(60, '*') + "\nError in file " + file
                           + ", line " + std::to_string(line) + ":\n\t" + x;
  OutFile << err_string << endl;
  throw std::runtime_error(err_string);
}
