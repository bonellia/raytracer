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
    friend parser::Vec3f operator+(const parser::Vec3f &lhs, const parser::Vec3f& rhs);
    friend parser::Vec3f operator-(const parser::Vec3f &lhs, const parser::Vec3f& rhs);
    unsigned char* InitializeImage(int width, int height);
    unsigned char * RenderScene(parser::Scene scene, int camera_no);
};

#endif // __raytracer_h__