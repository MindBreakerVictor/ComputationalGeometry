#ifndef CG_PCH_H
#define CG_PCH_H

#ifdef COMPUTATIONALGEOMETRY_EXPORTS
#define CGAPI __declspec(dllexport)
#else
#define CGAPI __declspec(dllimport)
#endif

#include <Windows.h>

#include <array>
#include <vector>


#endif

