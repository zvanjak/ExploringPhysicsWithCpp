#include "MMLBase.h"

void Demo_BaseUtils();
void Demo_Geometry();
void Demo_Geometry_2D();
void Demo_Geometry_3D();
void Demo_Matrix();
void Demo_MatrixNM();
void Demo_Vector();
//void Demo_VectorN();

int main()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                      EXAMPLE 1 - Basic objects                ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Demo_BaseUtils();
	Demo_Geometry();
	Demo_Geometry_2D();
	Demo_Geometry_3D();
	Demo_Matrix();
	Demo_MatrixNM();
	Demo_Vector();
	//Demo_VectorN();

  return 0;
}

