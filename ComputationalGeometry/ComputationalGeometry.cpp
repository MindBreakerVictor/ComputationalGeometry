#include "PCH.h"
#include "ComputationalGeometry.h"
#include <cmath>

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

Line2D::Line2D(const Point2D& X, const Point2D& Y)
{
	coefficients[0] = X.GetY() - Y.GetY();
	coefficients[1] = Y.GetX() - X.GetX();
	coefficients[2] = X.GetX() * Y.GetY() - Y.GetX() * X.GetY();
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

	if (Y == Z)
	{
		coefficients[0] = 0.0;
		coefficients[1] = 1.0;
		coefficients[2] = -1.0;
		return true;
	}

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

bool ComputationalGeometry::GetIntersection(const Segment2D& seg1, const Segment2D& seg2, Point2D** intersectionPoint, Segment2D** intersectionSegment)
{
	if (double delta = seg1.GetCoef()[0] * seg2.GetCoef()[1] - seg1.GetCoef()[1] * seg2.GetCoef()[0])
	{
		*intersectionPoint = new Point2D((-seg1.GetCoef()[2] * seg2.GetCoef()[1] + seg1.GetCoef()[1] * seg2.GetCoef()[2]) / delta, 
			(-seg1.GetCoef()[0] * seg2.GetCoef()[2] + seg1.GetCoef()[2] * seg2.GetCoef()[0]) / delta);

		Pair<double, double> interval = std::make_pair(min(seg1.GetA().GetY(), seg1.GetB().GetY()), max(seg1.GetA().GetY(), seg1.GetB().GetY()));

		if (seg1.GetA().GetX() == seg1.GetA().GetX() && interval.first <= (*intersectionPoint)->GetY() && interval.second >= (*intersectionPoint)->GetY())
		{
			interval = std::make_pair(min(seg2.GetA().GetX(), seg2.GetB().GetX()), max(seg2.GetA().GetX(), seg2.GetB().GetX()));

			if (interval.first <= (*intersectionPoint)->GetX() && interval.second >= (*intersectionPoint)->GetX())
				return true;

			delete *intersectionPoint;
			return false;
		}

		interval = std::make_pair(min(seg1.GetA().GetX(), seg1.GetB().GetX()), max(seg1.GetA().GetX(), seg1.GetB().GetX()));

		if (interval.first <= (*intersectionPoint)->GetX() && interval.second >= (*intersectionPoint)->GetX())
		{
			interval = std::make_pair(min(seg2.GetA().GetY(), seg2.GetB().GetY()), max(seg2.GetA().GetY(), seg2.GetB().GetY()));

			if (seg1.GetA().GetY() == seg2.GetB().GetY() && interval.first <= (*intersectionPoint)->GetY() && interval.second >= (*intersectionPoint)->GetY())
				return true;

			interval = std::make_pair(min(seg2.GetA().GetX(), seg2.GetB().GetX()), max(seg2.GetA().GetX(), seg2.GetB().GetX()));

			if (interval.first <= (*intersectionPoint)->GetX() && interval.second >= (*intersectionPoint)->GetY())
				return true;

			delete *intersectionPoint;
			return false;
		}

		delete *intersectionPoint;
		return false;
	}
	
	if (seg1.GetCoef()[0] * seg2.GetCoef()[2] - seg1.GetCoef()[2] * seg2.GetCoef()[0] == 0 &&
		seg1.GetCoef()[1] * seg2.GetCoef()[2] - seg1.GetCoef()[2] * seg2.GetCoef()[1] == 0)
	{
		if (seg1.GetA().GetX() < seg2.GetA().GetX() && seg1.GetA().GetX() < seg2.GetB().GetX() &&
			seg1.GetB().GetX() < seg2.GetA().GetX() && seg1.GetB().GetX() < seg2.GetB().GetX())
			return false;

		Point2D left = seg1.GetA().GetX() < seg2.GetA().GetX() ? seg2.GetA() : seg1.GetA();
		Point2D right = seg1.GetB().GetX() < seg2.GetB().GetX() ? seg1.GetB() : seg2.GetB();
		*intersectionSegment = new Segment2D(left, right);
		return true;
	}

	return false;
}

bool ComputationalGeometry::IsConvex(const Point2D& A, const Point2D& B, const Point2D& C, const Point2D& D)
{
	auto convexDeterminer = [](const Point2D& A, const Point2D& B, double x, double y)
	{
		return x * A.GetY() + A.GetX() * B.GetY() + B.GetX() * y - A.GetY() * B.GetX() - B.GetY() * x - y * A.GetX();
	};

	return convexDeterminer(A, B, C.GetX(), C.GetY()) * convexDeterminer(A, B, D.GetX(), D.GetY()) > 0.0 &&
		convexDeterminer(B, C, A.GetX(), A.GetY()) * convexDeterminer(B, C, D.GetX(), D.GetY()) > 0.0 &&
		convexDeterminer(C, D, A.GetX(), A.GetY()) * convexDeterminer(C, D, B.GetX(), B.GetY()) > 0.0 &&
		convexDeterminer(D, A, C.GetX(), C.GetY()) * convexDeterminer(D, A, B.GetX(), B.GetY()) > 0.0;
}

bool ComputationalGeometry::BelongsToConvexCoverage(const Point2D& A, const Point2D& B, const Point2D& C, const Point2D& D, const Point2D& M)
{
	auto area = [](const Point2D& A, const Point2D& B, const Point2D& C)
	{
		return abs(A.GetX() * B.GetY() + C.GetX() * A.GetY() + B.GetX() * C.GetY() -
			C.GetX() * B.GetY() - B.GetX() * A.GetY() - A.GetX() * C.GetY()) / 2.0;
	};

	return area(A, B, C) == area(M, A, B) + area(M, A, C) + area(M, B, C) ||
		area(A, D, C) == area(M, A, D) + area(M, A, C) + area(M, D, C) ||
		area(A, D, B) == area(M, A, D) + area(M, A, B) + area(M, D, B) ||
		area(C, D, B) == area(M, C, D) + area(M, C, B) + area(M, D, B);
}

bool ComputationalGeometry::IsDistrict(const Point2D& A, const Point2D& B, const Point2D& C, const Point2D& D)
{
	return GetDistance(A, B) + GetDistance(C, D) == GetDistance(A, D) + GetDistance(B, C);
}

double ComputationalGeometry::GetDistance(const Point2D& A, const Point2D& B)
{
	return sqrt(pow(B.GetX() - A.GetX(), 2) + pow(B.GetY() - A.GetY(), 2));
}

Position ComputationalGeometry::GetTCCRPosition(const Triangle2D& triangle, const Point2D& D)
{
	Array<Point2D, 3> trianglePoints = triangle.GetPoints();

	Pair<double, double> a(trianglePoints[0].GetX() - trianglePoints[1].GetX(), trianglePoints[0].GetY() - trianglePoints[1].GetY());
	Pair<double, double> b(trianglePoints[2].GetX() - trianglePoints[1].GetX(), trianglePoints[2].GetY() - trianglePoints[1].GetY());

	double cosB = (a.first * b.first + a.second * b.second) /
		(sqrt(pow(a.first, 2) + pow(a.second, 2)) * sqrt(pow(b.first, 2) + pow(b.second, 2)));

	a = std::make_pair(trianglePoints[2].GetX() - D.GetX(), trianglePoints[2].GetY() - D.GetY());
	b = std::make_pair(trianglePoints[0].GetX() - D.GetX(), trianglePoints[0].GetY() - D.GetY());

	double cosD = (a.first * b.first + a.second * b.second) /
		(sqrt(pow(a.first, 2) + pow(a.second, 2)) * sqrt(pow(b.first, 2) + pow(b.second, 2)));

	double angleSum = (acos(cosB) * 180.0 / PI) + (acos(cosD) * 180.0 / PI);

	if (angleSum == 180.0)
		return TCCRP_ONCIRCLE;
	else if (angleSum < 180.0)
		return TCCRP_OUTSIDE;

	return TCCRP_INSIDE;
}

Orientation ComputationalGeometry::GetOrientation(const Point2D& A, const Point2D& B, const Point2D& C)
{
	double value = (B.GetY() - A.GetY()) * (C.GetX() - B.GetX()) - (B.GetX() - A.GetX()) * (C.GetY() - B.GetY());

	if (!value)
		return ORIENTATION_COLLINEAR;

	return value > 0.0 ? ORIENTATION_CLOCKWISE : ORIENTATION_COUNTERCLOCKWISE;
}

Vector<Point2D> ComputationalGeometry::GetConvexHullBorder(const Vector<Point2D>& convexHullPoints)
{
	if (convexHullPoints.size() < 3)
		return convexHullPoints;

	size_t bottomLeftPointIndex = 0;

	// Get the bottom left point.
	for (size_t i = 1; i < convexHullPoints.size(); ++i)
		if (convexHullPoints[bottomLeftPointIndex].GetX() == convexHullPoints[i].GetX() &&
			convexHullPoints[bottomLeftPointIndex].GetY() > convexHullPoints[i].GetY() ||
			convexHullPoints[bottomLeftPointIndex].GetX() > convexHullPoints[i].GetX())
				bottomLeftPointIndex = i;

	Vector<Point2D> convexHullBorder;

	// Start from the bottom left point, moving counterclockwise until we reach the start point.
	int currentPointIndex = bottomLeftPointIndex;

	do 
	{
		convexHullBorder.push_back(convexHullPoints[currentPointIndex]);

		int mostCounterclockwisePointIndex = (currentPointIndex + 1) % convexHullPoints.size();
		
		for (size_t i = 0; i < convexHullPoints.size(); ++i)
			if (GetOrientation(convexHullPoints[currentPointIndex], convexHullPoints[i],
				convexHullPoints[mostCounterclockwisePointIndex]) == ORIENTATION_COUNTERCLOCKWISE)
				mostCounterclockwisePointIndex = i;

		currentPointIndex = mostCounterclockwisePointIndex;
	} while (currentPointIndex != bottomLeftPointIndex);

	return convexHullBorder;
}

