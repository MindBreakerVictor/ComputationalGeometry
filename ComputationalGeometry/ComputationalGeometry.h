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

	/* Constants */
	const Point2D Origin2D;
	const Point3D Origin3D;
}

#endif
