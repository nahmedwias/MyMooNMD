// =======================================================================
// %W% %G%
//
// Class:       TFEDesc3D
// Purpose:     store a finite element descriptor for a 3D element
//
// Author:      Gunar Matthies  23.07.98
//
// =======================================================================

#include <FEDesc3D.h>
#include <string.h>
#include <iostream>

/** constructor, setting all data */
TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof,
                     int n_outerdof, int *outerdof)
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = n_outerdof;
  OuterDOF = outerdof;
  
#ifdef _MPI    
  N_EdgeDOF = 0;
  EdgeDOF = nullptr;
  N_VertDOF = 0;
  VertDOF = nullptr;
  EdgeVertData_Filled = 0;
#endif
  
}


#ifdef _MPI 
TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof,
                     int n_edgeDOF, int **edgeDOF, int n_vertDOF, int *vertDOF )
 : TFEDesc3D(description, n_dof, n_jointdof, jointdof, n_innerdof, innerdof,
             0, nullptr, n_edgeDOF, edgeDOF, n_vertDOF, vertDOF)
{
}


TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof,
                     int n_outerdof, int *outerdof,
                     int n_edgeDOF, int **edgeDOF, int n_vertDOF, int *vertDOF )
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = n_outerdof;
  OuterDOF = outerdof;

  N_EdgeDOF = n_edgeDOF;
  EdgeDOF = edgeDOF;
  N_VertDOF = n_vertDOF;
  VertDOF = vertDOF; 
  EdgeVertData_Filled = 1;  
}

#endif // _MPI

/** return face on which the i-th local degree of freedom is   
If i is not a dof on a face, return -1

If i is a dof on two faces (e.g. on a vertex), one of these two faces is 
returned. Don't use this function in this case.

*/
int TFEDesc3D::GetJointOfThisDOF(int localDOF) const
{
  bool is_DOF_on_edge = false;
 
  for (int i = 0; i < N_OuterDOF; i++)
  {
    if(OuterDOF[i]==localDOF)
    {
      is_DOF_on_edge=true;
      break;
    }
  }
  if(!is_DOF_on_edge)
    return -1;
  //else // continue to find the edge
  int i = 0;
  while (true)
  {
    // this must terminate, since we already know localDOF is a dof on an edge
    for (int j = 0; j < N_JointDOF; j++)
    {
      if(JointDOF[i][j] == localDOF)
        return i;
    }
    i++;
  }
}
