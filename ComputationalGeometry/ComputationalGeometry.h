#ifndef COMPUTATIONAL_GEOMETRY_H
#define COMPUTATIONAL_GEOMETRY_H

#include "PCH.h"

namespace ComputationalGeometry
{
	template <class Ty1, class Ty2>
	using Pair = std::pair<Ty1, Ty2>;

	template <class Ty, size_t size>
	using Array = std::array<Ty, size>;

	template <class Ty, class Alloc = std::allocator<Ty>>
	using Vector = std::vector<Ty, Alloc>;

	class CGAPI Point2D
	{
		public:
			Point2D() : _x(), _y() { }
			Point2D(double x, double y) : _x(x), _y(y) { }
			Point2D(const Point2D& source) : _x(source._x), _y(source._y) { }

			double GetX() const { return _x; }
			double GetY() const { return _y; }

			void SetX(double x) { _x = x; }
			void SetY(double y) { _y = y; }

			Point2D& operator=(const Point2D& source) = default;

			Point2D operator-(const Point2D& right) const { return Point2D(_x - right._x, _y - right._y); }
			Point2D operator+(const Point2D& right) const { return Point2D(_x + right._x, _y + right._y); }

			Point2D operator*(double scalar) const { return Point2D(_x * scalar, _y * scalar); }
			Point2D operator/(double scalar) const { return Point2D(_x / scalar, _y / scalar); }

			Point2D& operator+=(const Point2D& right);
			Point2D& operator-=(const Point2D& right);
			Point2D& operator*=(double scalar);
			Point2D& operator/=(double scalar);

			bool operator==(const Point2D& right) const { return _x == right._x && _y == right._y; }
			bool operator!=(const Point2D& right) const { return !(*this == right); }

		private:
			double _x;
			double _y;
	};

	class CGAPI Point3D
	{
		public:
			Point3D() : _x(), _y() { }
			Point3D(double x, double y, double z) : _x(x), _y(y), _z(z) { }
			Point3D(const Point3D& source) : _x(source._x), _y(source._y), _z(source._z) { }

			double GetX() const { return _x; }
			double GetY() const { return _y; }
			double GetZ() const { return _z; }

			void SetX(double x) { _x = x; }
			void SetY(double y) { _y = y; }
			void SetZ(double z) { _z = z; }

			Point3D& operator=(const Point3D& source) = default;

			Point3D operator-(const Point3D& right) const { return Point3D(_x - right._x, _y - right._y, _z - right._z); }
			Point3D operator+(const Point3D& right) const { return Point3D(_x + right._x, _y + right._y, _z + right._z); }

			Point3D operator*(double scalar) const { return Point3D(_x * scalar, _y * scalar, _z * scalar); }
			Point3D operator/(double scalar) const { return Point3D(_x / scalar, _y / scalar, _z / scalar); }

			Point3D& operator+=(const Point3D& right);
			Point3D& operator-=(const Point3D& right);
			Point3D& operator*=(double scalar);
			Point3D& operator/=(double scalar);

			bool operator==(const Point3D& right) const { return _x == right._x && _y == right._y && _z == right._z; }
			bool operator!=(const Point3D& right) const { return !(*this == right); }

		private:
			double _x;
			double _y;
			double _z;
	};

	class CGAPI Line2D
	{
		public:
			Line2D(const Point2D& X, const Point2D& Y);

			const Array<double, 3>& GetCoef() const { return coefficients; }

		protected:
			Array<double, 3> coefficients;
	};

	class CGAPI Segment2D : public Line2D
	{
		public:
			Segment2D(const Point2D& X, const Point2D& Y) : A(X), B(Y), Line2D(X, Y) { }

			const Point2D& GetA() const { return A; }
			const Point2D& GetB() const { return B; }

		private:
			Point2D A;
			Point2D B;
	};

	bool CGAPI Collinear(const Point2D& X, const Point2D& Y, const Point2D& Z);
	bool CGAPI Collinear(const Point3D& X, const Point3D& Y, const Point3D& Z);

	// If the points passed as parameters are collinear 
	// than the vectors XY and YZ can be written as XY = r * YZ
	// The function returns false if they aren't collinear and the value of ratio remains unchanged.
	// Otherwise it returns true. If the ratio isn't equal to -1 than there is only one solution.
	// If the ratio is -1 it means there are an infinite number of solutions.
	bool CGAPI GetRatio(const Point3D& X, const Point3D& Y, const Point3D& Z, double& ratio);

	// If the points passed as parameters are collinear the function returns true.
	// The barycentric combination coefficients are inserted in the array passed in coefficients parameter 
	// only if there is an unique solution, a call to GetRation returns true and parameter ratio isn't equal to -1).
	bool CGAPI GetBarycentricCombination(const Point3D& X, const Point3D& Y, const Point3D& Z, Array<double, 3>& coefficients);

	/*
		Pass the address of a Point2D pointer(the pointer should be initialized with nullptr) in intersectionPoint parameter.
		Pass the address of a Segment2D pointer(the pointer should be initialized with nullptr) in intersectionSegment parameter.
		If the function returns false the two segments passed in arguments seg1 & seg2 do not intersect.
		If the function returns true there are two possible cases: 
			*intersectionPoint points to a valid Point2D which is the intersection point of the two segments.
			*intersectionSegment points to a valid Segment2D which is the interval of points where the two segments intersect.
	*/
	bool CGAPI GetIntersection(const Segment2D& seg1, const Segment2D& seg2, Point2D** intersectionPoint, Segment2D** intersectionSegment);

	// Check if the quadrilateral formed by the four points passed is convex.
	bool CGAPI IsConvex(const Point2D& A, const Point2D& B, const Point2D& C, const Point2D& D);

	// Check if the point M belongs to the convex coverage of points A, B, C, D.
	bool CGAPI BelongsToConvexCoverage(const Point2D& A, const Point2D& B, const Point2D& C, const Point2D& D, const Point2D& M);

	// Check if the quadrilateral formed by the four points passed as arguments is district.
	bool CGAPI IsDistrict(const Point2D& A, const Point2D& B, const Point2D& C, const Point2D& D);

	// Returns the distance between two 2D points.
	double CGAPI GetDistance(const Point2D& A, const Point2D& B);

	/* Constants */
	const Point2D Origin2D;
	const Point3D Origin3D;
}

#endif

