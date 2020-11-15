#ifndef __rayrtacer_h__
#define __raytracer_h__

#include "parser.h"
#include "definitions.h"
#include <cmath>

struct Ray {
    Vec3f origin;
    Vec3f direction;
};

struct Touch {
    float t;
    parser::Vec3f position;
    parser::Vec3f normal;
    int material_id;

};

static constexpr Touch MISS = {-1, {0, 0, 0}, {0, 0, 0}, -1};

class RayTracer {
public:

    Scene scene;

    // Geometry related methods
    Vec3f Cross(const Vec3f &lhs, const Vec3f &rhs);

    float Dot(const Vec3f &lhs, const Vec3f &rhs);

    // Multiplies a vector with the given float number on left-hand-side.
    Vec3f Scale(const float &lhs, const Vec3f &rhs);

    // Adds two vectors.
    Vec3f Add(const parser::Vec3f &lhs, const Vec3f &rhs);

    // Subtracts two vectors.
    Vec3f Subtract(const parser::Vec3f &lhs, const Vec3f &rhs);

    // Negates a vector.
    Vec3f Negate(const parser::Vec3f &rhs);

    // Calculates the length of a vector.
    float Length(parser::Vec3f vector);

    // Normalizes a vector.
    Vec3f Normalize(Vec3f vector);

    // Generates a ray from camera to near plane.
    Ray GenerateEyeRay(int x, int y, Camera cam);

    // Tests a ray intersection with the given sphere, returns Touch information.
    Touch IntersectSphere(Ray ray, Sphere sphere);

    // Fills the initial image with red colors.
    static unsigned char *InitializeImage(int width, int height);

    // Renders the scene.
    static unsigned char *RenderScene(Scene scene, int camera_no);
};

#endif // __raytracer_h__