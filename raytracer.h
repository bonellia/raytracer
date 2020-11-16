#ifndef __rayrtacer_h__
#define __raytracer_h__

#include "parser.h"
#include "definitions.h"
#include <cmath>

struct Ray {
    Vec3f origin;
    Vec3f direction;
};
/*
    For debugging purposes.
    Sphere, Triangle, Mesh, None
*/
enum Contact {
    BUM, BOOB, BODY, NAH
};

struct TouchAttempt {
    float t;
    parser::Vec3f position;
    parser::Vec3f normal;
    int material_id;
    Contact contact;
    int touched_object_no;
};

static constexpr TouchAttempt MISS = {-1, {0, 0, 0}, {0, 0, 0}, -1,  NAH, -1};

class RayTracer {
public:

    Scene scene;
    std::vector<TouchAttempt> touch_list;

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
    Vec3f Add(const Vec3f &lhs, const Vec3f &rhs);

    // Subtracts two vectors.
    Vec3f Subtract(const Vec3f &lhs, const Vec3f &rhs);

    // Negates a vector.
    Vec3f Negate(const Vec3f &vector);

    // Calculates the length of a vector.
    float Length(const Vec3f &vector);

    // Normalizes a vector.
    Vec3f Normalize(const Vec3f &vector);

    float Determinant(const Vec3f &a, const Vec3f &b, const Vec3f &c);

    // Generates a ray from camera to near plane.
    Ray GenerateEyeRay(const int pixel_row, const int pixel_column, const Camera &cam);

    // Tests a ray intersection with the given sphere, returns TouchAttempt information.
    TouchAttempt SphereIntersectionTest(const Ray &ray, const Sphere &sphere);

    // Tests a ray intersection with the given triangle, returns TouchAttempt information.
    TouchAttempt TriangleIntersectionTest(const Ray &ray, const Vec3f &a, const Vec3f &b, const Vec3f &c, const int material_id);

    // Tests ray intersections for the triangles of a mesh and updates touch list.
    TouchAttempt MeshIntersectionTest(const Ray &ray, const Mesh &mesh);

    // Tests ray intersections for all objects within the scene and returns the closest touching point information.
    TouchAttempt FindClosestContact(const Ray &ray);

    // Fills the initial image with red colors.
    static unsigned char *InitializeImage(const int width, const int height);

    // Calculates the color value for a pixel considering all scene parameters.
    Vec3f CalculatePixelColor(const Ray &ray, int depth);

    // Renders the scene for the given camera.
    unsigned char *RenderScene(const Camera &camera, const int width, const int height);
};

#endif // __raytracer_h__