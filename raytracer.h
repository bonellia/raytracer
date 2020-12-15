#ifndef __raytracer_h__
#define __raytracer_h__

#include "parser.h"
#include "definitions.h"
#include <cmath>

typedef unsigned char RGB[3];

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

static constexpr TouchAttempt MISS = {-1, {0, 0, 0}, {0, 0, 0}, -1, NAH, -1};

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
    Ray GenerateEyeRay(int pixel_row, int pixel_column, const Camera &cam);

    // Tests a ray intersection with the given sphere, returns TouchAttempt information.
    TouchAttempt SphereIntersectionTest(const Ray &ray, const Sphere &sphere);

    // Tests a ray intersection with the given triangle, returns TouchAttempt information.
    TouchAttempt
    TriangleIntersectionTest(const Ray &ray, const Vec3f &a, const Vec3f &b, const Vec3f &c, int material_id);

    // Tests ray intersections for the triangles of a mesh and updates touch list.
    TouchAttempt MeshIntersectionTest(const Ray &ray, const Mesh &mesh);

    // Tests ray intersections for all objects within the scene and returns the closest touching point information.
    TouchAttempt FindClosestContact(const Ray &ray);

    /*!
     * Fills the initial image with red colors. Intended for all image.
     * @param width Column count as pixels.
     * @param height Row count as pixels.
     * @return
     */
    static unsigned char *InitializeImage(int width, int height);

    /*!
     * Calculates the color value for a pixel considering all scene parameters.
     * @param ray Ray that intersects the current pixel value on near plane.
     * @param touch_attempt Information for the ray surface intersection test.
     * @param cam Camera that used for color calculation. Needed for external rays.
     * @param depth Maximum number of recursive calls.
     * @return The vector with color data.
     */
    Vec3f CalculatePixelColor(const Ray &ray, const Camera &cam, int depth);

    /*!
     * * Maps RGB values to image array.
     * @param image The image.
     * @param row Row number of the pixel.
     * @param column Column number of the pixel.
     * @param color RGB value as a triplet.
     * @param width Total width of the image.
     * @param color A triplet with red, green, blue values.
     */
    void SetImagePixelRGB(unsigned char *&image, int row, int column, int width, const RGB color);
    /*!
     * Utilizes parallelization by rending a vertical bar fraction of the desired image.
     * @param cam Camera that used for rendering.
     * @param image The image.
     * @param width_from Left bound index of the image part.
     * @param width_to Right bound index of the image part.
     * @param height Row count to render.
     */
    void RenderBar(const Camera &cam, unsigned char *&image, int width_from, int width_to, int height);

    /*!
     * Renders the scene for the given camera.
     * @param cam Camera that used for rendering.
     * @param width Column count to render.
     * @param height Row count to render.
     * @return One dimensional array that contains the color values (e.g., 255255255...)
     */
    unsigned char *RenderScene(const Camera &cam, int width, int height);
};

#endif // __raytracer_h__