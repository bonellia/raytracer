#ifndef __rayrtacer_h__
#define __raytracer_h__

#include "parser.h"
#include <cmath>
class RayTracer
{
public:
    parser::Scene scene;
    // Geometry related methods
    parser::Vec3f Cross(const parser::Vec3f& lhs, const parser::Vec3f& rhs);
    float Dot(const parser::Vec3f& lhs, const parser::Vec3f& rhs);
    // Multiplies a vector with the given float number on left-hand-side.
    friend parser::Vec3f operator*(const float &lhs, const parser::Vec3f &rhs);
    friend parser::Vec3f operator+(const parser::Vec3f &lhs, const parser::Vec3f& rhs);
    friend parser::Vec3f operator-(const parser::Vec3f &lhs, const parser::Vec3f& rhs);
    float Length(parser::Vec3f vector);
    parser::Vec3f Normalize(parser::Vec3f vector);
    parser::Ray GenerateEyeRay(int x, int y, parser::Camera cam);
    unsigned char* InitializeImage(int width, int height);
    unsigned char * RenderScene(parser::Scene scene, int camera_no);
};

#endif // __raytracer_h__