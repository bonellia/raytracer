#include "util.h"
#include <iostream>
/*
Overloads << operator for Vec3f object, allowing it to be more human readable.
*/
std::ostream& operator<<(std::ostream& os, const parser::Vec3f& vec)
{
    os << '(' << vec.x << ", " << vec.y << ", " << vec.z << ')';
    return os;
}
/*
Overloads << operator for Vec3i object, allowing it to be more human readable.
*/
std::ostream& operator<<(std::ostream& os, const parser::Vec3i& vec)
{
    os << '(' << vec.x << ", " << vec.y << ", " << vec.z << ')';
    return os;
}
/*
Overloads << operator for Vec4f object, allowing it to be more human readable.
*/
std::ostream& operator<<(std::ostream& os, const parser::Vec4f& plane)
{
    os << '(' << plane.x << ", " << plane.y << ", " << plane.z << ", " << plane.w << ')';
    return os;
}
/*
Overloads << operator for Camera object, allowing it to be more human readable.
*/
std::ostream& operator<<(std::ostream& os, const parser::Camera& cam)
{
    os << "Image name: " << cam.image_name << '\n'
       << "H:"<< cam.image_height << ", W: "<< cam.image_width << '\n'
       << "Near plane: " << cam.near_plane << '\n'
       << "Near distance:" << cam.near_distance << '\n'
       << "Position: " << cam.position << '\n'
       << "Up vector: " << cam.up << '\n'
       << "Gaze vector: " << cam.gaze << '\n';
    return os;
}
/*
Overloads << operator for PointLight object, allowing it to be more human readable.
*/
std::ostream& operator<<(std::ostream& os, const parser::PointLight& pl)
{
    os << "Intensity: " << pl.intensity << '\n'
       << "Position:"<< pl.position << '\n';
    return os;
}
/*
Overloads << operator for Material object, allowing it to be more human readable.
*/
std::ostream& operator<<(std::ostream& os, const parser::Material& mat)
{
    os << "Ambient: " << mat.ambient << '\n'
       << "Diffuse: " << mat.diffuse << '\n'
       << "Specular: " << mat.specular << '\n'
       << "Mirror: " << mat.mirror << '\n'
       << "Phong exponent: " << mat.phong_exponent << '\n';
    return os;
}
/*
Overloads << operator for Sphere object, allowing it to be more human readable.
*/
std::ostream& operator<<(std::ostream& os, const parser::Sphere& sph)
{
    os << "Material ID: " << sph.material_id << '\n'
       << "Center vertex ID: " << sph.center_vertex_id << '\n'
       << "Radius: " << sph.radius << '\n';
    return os;
}
/*
Overloads << operator for Scene object, allowing it to be more human readable.
*/
std::ostream& operator<<(std::ostream& os, const parser::Scene& scene)
{
    int counter = 0;
    os
            << "------------------------------\n"
            << "       SCENE PROPERTIES       \n"
            << "------------------------------\n"
            << "Background color: " << scene.background_color << '\n'
            << "Shadow ray epsilon: " << scene.shadow_ray_epsilon << '\n'
            << "Maximum recursion depth: " << scene.max_recursion_depth << '\n'
            << "Ambient light: " << scene.ambient_light << '\n'
            << "Vertex count: " << scene.vertex_data.size() << '\n'
            << "Mesh count: " << scene.meshes.size() << '\n'
            << "Triangle count: " << scene.triangles.size() << '\n'
            << "------------------------------\n"
            << "          Camera(s)           \n"
            << "------------------------------\n";
    for (parser::Camera c : scene.cameras)
    {
        os << "Camera " << counter++ << ":\n---------\n" << c;
    }
    counter = 1;
    os
            << "------------------------------\n"
            << "        Point light(s)        \n"
            << "------------------------------\n";
    for (parser::PointLight pl : scene.point_lights)
    {
        os << "Point Light " << counter++ << ":\n--------------\n" << pl;
    }
    counter = 1;
    os
            << "------------------------------\n"
            << "         Material(s)          \n"
            << "------------------------------\n";
    for (parser::Material mat : scene.materials)
    {
        os << "Material " << counter++ << ":\n-----------\n" << mat;
    }
    counter = 1;
    os
            << "------------------------------\n"
            << "           Spheres            \n"
            << "------------------------------\n";
    for (parser::Sphere sph : scene.spheres)
    {
        os << "Sphere " << counter++ << ":\n---------\n" << sph;
    }
    return os;

}
/*
Prints the scene configuration in a human readable format.
*/
void util::Util::PrintSceneDetails(parser::Scene scene)
{
    std::cout << scene << std::endl;
}