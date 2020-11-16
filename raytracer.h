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
    int touched_object_no;
};

static constexpr Touch MISS = {-1, {0, 0, 0}, {0, 0, 0}, -1, -1};

class RayTracer {
public:

    Scene scene;
    std::vector<Touch> touch_list;

    // Geometry related methods
    Vec3f Cross(const Vec3f &lhs, const Vec3f &rhs);

    float Dot(const Vec3f &lhs, const Vec3f &rhs);

    // Multiplies rhs vector with the given float number on left-hand-side.
    Vec3f Scale(const float &lhs, const Vec3f &rhs);

    // Multiplies rhs vector with the components of the lhs vector.
    Vec3f VectorScale(const Vec3f &lhs, const Vec3f &rhs);

    // Divides lhs vector with the given float number on right-hand-side.
    Vec3f Circumsize(const Vec3f &lhs, const float &rhs);

    // Adds two vectors.
    Vec3f Add(const parser::Vec3f &lhs, const Vec3f &rhs);

    // Subtracts two vectors.
    Vec3f Subtract(const parser::Vec3f &lhs, const Vec3f &rhs);

    // Negates a vector.
    Vec3f Negate(const parser::Vec3f &vector);

    // Calculates the length of a vector.
    float Length(parser::Vec3f vector);

    // Normalizes a vector.
    Vec3f Normalize(Vec3f vector);

    float Determinant(Vec3f a, Vec3f b, Vec3f c);

    // Generates a ray from camera to near plane.
    Ray GenerateEyeRay(int pixel_row, int pixel_column, Camera cam);

    // Tests a ray intersection with the given sphere, returns Touch information.
    Touch SphereIntersectionTest(Ray ray, Sphere sphere);

    // Tests a ray intersection with the given triangle, returns Touch information.
    Touch TriangleIntersectionTest(Ray ray, Vec3f &a, Vec3f &b, Vec3f &c, int material_id);

    // Tests ray intersections for the triangles of a mesh and updates touch list.
    Touch MeshIntersectionTest(Ray ray, Mesh mesh);

    // Tests ray intersections for all objects within the scene and returns the closest touching point information.
    Touch FindClosestTouch(Ray ray);

    // Fills the initial image with red colors.
    static unsigned char *InitializeImage(int width, int height);

    // Calculates the color value for a pixel considering all scene parameters.
    Vec3f CalculatePixelColor(Ray ray, int depth);

    // Renders the scene for the given camera.
    unsigned char *RenderScene(Camera camera, int width, int height);
};

#endif // __raytracer_h__