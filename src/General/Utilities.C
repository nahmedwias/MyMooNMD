#include <Utilities.h>
#include <chrono>
#include <ctime>
#include <unistd.h>

std::string utilities::get_host_name()
{
  char buf[80];
  gethostname(buf,80);
  return std::string(buf);
}


std::string utilities::get_date_and_time()
{
  auto time = std::chrono::system_clock::now();
  std::time_t t = std::chrono::system_clock::to_time_t(time);
  return std::ctime(&t);
}
