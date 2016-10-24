#include "PCH.h"
#include "ComputationalGeometry.h"

using namespace ComputationalGeometry;

Point2D& Point2D::operator+=(const Point2D& right)
{
	_x += right._x;
	_y += right._y;
	return *this;
}

Point2D& Point2D::operator-=(const Point2D& right)
{
	_x -= right._x;
	_y -= right._y;
	return *this;
}

Point2D& Point2D::operator*=(double scalar)
{
	_x *= scalar;
	_y *= scalar;
	return *this;
}

Point2D& Point2D::operator/=(double scalar)
{
	_x /= scalar;
	_y /= scalar;
	return *this;
}

Point3D& Point3D::operator+=(const Point3D& right)
{
	_x += right._x;
	_y += right._y;
	_z += right._z;
	return *this;
}

Point3D& Point3D::operator-=(const Point3D& right)
{
	_x -= right._x;
	_y -= right._y;
	_z -= right._z;
	return *this;
}

Point3D& Point3D::operator*=(double scalar)
{
	_x *= scalar;
	_y *= scalar;
	_z *= scalar;
	return *this;
}

Point3D& Point3D::operator/=(double scalar)
{
	_x /= scalar;
	_y /= scalar;
	_z /= scalar;
	return *this;
}

bool ComputationalGeometry::Collinear(const Point2D& X, const Point2D& Y, const Point2D& Z)
{
	return ((Y.GetX() - X.GetX()) * (Z.GetY() - X.GetY()) - (Z.GetX() - X.GetX()) * (Y.GetY() - X.GetY())) == 0;
}

bool ComputationalGeometry::Collinear(const Point3D& X, const Point3D& Y, const Point3D& Z)
{
	return (((Y.GetY() - X.GetY()) * (Z.GetZ() - X.GetZ()) - (Z.GetY() - X.GetY()) * (Y.GetZ() - X.GetZ())) == 0.0 &&
		((Z.GetX() - X.GetX()) * (Y.GetZ() - X.GetZ()) - (Y.GetX() - X.GetX()) * (Z.GetZ() - X.GetZ())) == 0.0 &&
		((Y.GetX() - X.GetX()) * (Z.GetY() - X.GetY()) - (Z.GetX() - X.GetX()) * (Y.GetY() - X.GetY())) == 0.0);
}

bool ComputationalGeometry::GetRatio(const Point3D& X, const Point3D& Y, const Point3D& Z, double& ratio)
{
	Point3D XY(Y.GetX() - X.GetX(), Y.GetY() - X.GetY(), Y.GetZ() - X.GetZ());
	Point3D YZ(Z.GetX() - Y.GetX(), Z.GetY() - Y.GetY(), Z.GetZ() - Y.GetZ());

	if (XY == Origin3D && YZ == Origin3D)
	{
		ratio = -1.0;
		return true;
	}

	if (YZ == Origin3D)
		return false;

	if (XY == Origin3D)
	{
		ratio = 0.0;
		return true;
	}

	if (YZ.GetX() == 0.0)
	{
		if (XY.GetX() != 0.0)
			return false;

		if (YZ.GetY() == 0.0)
		{
			if (XY.GetY() != 0.0)
				return false;

			ratio = XY.GetZ() / YZ.GetZ();
			return true;
		}

		if (YZ.GetZ() == 0.0)
		{
			if (XY.GetZ() != 0.0)
				return false;

			ratio = XY.GetY() / YZ.GetY();
			return true;
		}

		double r = XY.GetY() / YZ.GetY();

		if (r == XY.GetZ() / YZ.GetZ())
		{
			ratio = r;
			return true;
		}

		return false;
	}

	if (YZ.GetY() == 0.0)
	{
		if (XY.GetY() != 0.0)
			return false;

		if (YZ.GetZ() == 0.0)
		{
			if (XY.GetZ() != 0.0)
				return false;

			ratio = XY.GetX() / YZ.GetX();
			return true;
		}

		double r = XY.GetX() / YZ.GetX();

		if (r == XY.GetZ() / YZ.GetZ())
		{
			ratio = r;
			return true;
		}

		return false;
	}

	double r = XY.GetX() / YZ.GetX();

	if (YZ.GetZ() == 0.0 && r == XY.GetY() / YZ.GetY())
	{
		ratio = r;
		return true;
	}

	if (r == XY.GetY() / YZ.GetY() && r == XY.GetZ() / YZ.GetZ())
	{
		ratio = r;
		return true;
	}

	return false;
}

bool ComputationalGeometry::GetBarycentricCombination(const Point3D& X, const Point3D& Y, const Point3D& Z, Array<double, 3>& coefficients)
{
	double ratio;

	if (!GetRatio(X, Y, Z, ratio))
		return false;

	if (ratio != -1.0)
	{
		coefficients[0] = 1.0;
		coefficients[1] = 1.0 + ratio;
		coefficients[2] = -ratio;
	}

	return true;
}

