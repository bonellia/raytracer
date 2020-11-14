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
    parser::Vec3f Multiply(const float &lhs, const parser::Vec3f &rhs);
    // Adds two vectors.
    parser::Vec3f Add(const parser::Vec3f &lhs, const parser::Vec3f& rhs);
    // Subtracts two vectors.
    parser::Vec3f Subtract(const parser::Vec3f &lhs, const parser::Vec3f& rhs);
    // Negates a vector.
    parser::Vec3f Negate(const parser::Vec3f& rhs);
    // Calculates the length of a vector.
    float Length(parser::Vec3f vector);
    // Normalizes a vector.
    parser::Vec3f Normalize(parser::Vec3f vector);
    // Generates a ray from camera to near plane.
    parser::Ray GenerateEyeRay(int x, int y, parser::Camera cam);
    // Fills the initial image with red colors.
    unsigned char* InitializeImage(int width, int height);
    // Renders the scene.
    unsigned char * RenderScene(parser::Scene scene, int camera_no);
};

#endif // __raytracer_h__