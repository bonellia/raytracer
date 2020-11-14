#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "raytracer.h"
#include "util.h"

typedef unsigned char RGB[3];

parser::Vec3f Cross(const parser::Vec3f &lhs, const parser::Vec3f &rhs)
{
    return {
            lhs.y * rhs.z - rhs.y * lhs.z,
            rhs.x * lhs.z - lhs.x * rhs.z,
            lhs.x * rhs.y - rhs.x * lhs.y};
}
float Dot(const parser::Vec3f &lhs, const parser::Vec3f &rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}
parser::Vec3f operator+(const parser::Vec3f &lhs, const parser::Vec3f &rhs)
{
    return {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z};
}
parser::Vec3f operator-(const parser::Vec3f &lhs, const parser::Vec3f &rhs)
{
    return {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z};
}

unsigned char* RayTracer::InitializeImage(int width, int height)
{
    unsigned char *image = new unsigned char[width * height * 3];
    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            image[i++] = 255;
            image[i++] = 0;
            image[i++] = 0;
        }
    }
    return image;
}

unsigned char * RayTracer::RenderScene(parser::Scene scene, int camera_no)
{
    int width = scene.cameras.at(camera_no).image_width;
    int height = scene.cameras.at(camera_no).image_height;
    unsigned char* raw_image = InitializeImage(width, height);
    // Step 1: Initialize the image.
    return raw_image;
}

int main(int argc, char *argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;
    util::Util util;
    RayTracer ray_tracer;
    scene.loadFromXml(argv[1]);
    ray_tracer.scene = scene;
    util.PrintSceneDetails(scene);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.
    /*
    const RGB BAR_COLOR[8] =
            {
                    {255, 255, 255}, // 100% White
                    {255, 255, 0},   // Yellow
                    {0, 255, 255},   // Cyan
                    {0, 255, 0},     // Green
                    {255, 0, 255},   // Magenta
                    {255, 0, 0},     // Red
                    {0, 0, 255},     // Blue
                    {0, 0, 0},       // Black
            };

    int width = 640, height = 480;
    int columnWidth = width / 8;

    unsigned char *image = new unsigned char[width * height * 3];

    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }
    */
    int width = scene.cameras.at(0).image_width;
    int height = scene.cameras.at(0).image_height;
    unsigned char *image = ray_tracer.RenderScene(scene, 0);
    write_ppm("test.ppm", image, width, height);
}
