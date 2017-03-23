#ifndef read_write_metis
#define read_write_metis

#include <mpi.h>
#include <Domain.h>
#include <metis.h>

class ReadWriteMetis
{

public:

	ReadWriteMetis() = default;
	~ReadWriteMetis() = default;

	void readFile(TDomain *Domain,int size, int N_Cells, idx_t *Cell_Rank);
	void writeFile(TDomain *Domain,int size, int N_Cells, idx_t *Cell_Rank);
};
#endif
