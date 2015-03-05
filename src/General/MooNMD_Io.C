#include <Database.h>
#include <MooNMD_Io.h>

std::ofstream OutFile;
std::ofstream LogFile;

void OpenFiles()
{
  OutFile.open(TDatabase::ParamDB->OUTFILE);
  LogFile.open(TDatabase::ParamDB->LOGFILE);
}

void CloseFiles()
{
  OutFile.close();
  LogFile.close();
}


