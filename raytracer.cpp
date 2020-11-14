#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "raytracer.h"
#include "util.h"

typedef unsigned char RGB[3];

parser::Vec3f RayTracer::Cross(const parser::Vec3f &lhs, const parser::Vec3f &rhs) {
    return {
            lhs.y * rhs.z - rhs.y * lhs.z,
            rhs.x * lhs.z - lhs.x * rhs.z,
            lhs.x * rhs.y - rhs.x * lhs.y};
}

float RayTracer::Dot(const parser::Vec3f &lhs, const parser::Vec3f &rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

parser::Vec3f RayTracer::Scale(const float &lhs, const parser::Vec3f &rhs) {
    parser::Vec3f vector;
    vector.x = lhs * rhs.x;
    vector.y = lhs * rhs.y;
    vector.z = lhs * rhs.z;
    return vector;
}

parser::Vec3f RayTracer::Add(const parser::Vec3f &lhs, const parser::Vec3f &rhs) {
    return {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z};
}

parser::Vec3f RayTracer::Subtract(const parser::Vec3f &lhs, const parser::Vec3f &rhs) {
    return {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z};
}

parser::Vec3f RayTracer::Negate(const parser::Vec3f &rhs) {
    return {
            -rhs.x,
            -rhs.y,
            -rhs.z};
}


float RayTracer::Length(parser::Vec3f vector) {
    return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

parser::Vec3f RayTracer::Normalize(parser::Vec3f vector) {
    parser::Vec3f normalized_vector;
    float vector_length = Length(vector);
    normalized_vector.x = vector.x / vector_length;
    normalized_vector.y = vector.y / vector_length;
    normalized_vector.z = vector.z / vector_length;
    return normalized_vector;
}

parser::Ray RayTracer::GenerateEyeRay(int x, int y, parser::Camera cam) {
    parser::Ray ray;
    parser::Vec3f su, sv, s;
    ray.origin = cam.position;
    float pixel_width = cam.near_plane.y - cam.near_plane.x / (float) cam.image_width;
    float pixel_height = (cam.near_plane.w - cam.near_plane.z) / cam.image_height;
    // TODO: continue with su & sv.
    parser::Vec3f cam_u = Cross(cam.gaze, cam.up);
    su = Scale(cam.near_plane.x + (x + 0.5) * (pixel_width), cam_u);
    sv = Scale(cam.near_plane.z + (x + 0.5) * (pixel_height), cam.up);
    s = Add(su, sv);

    ray.direction = Add(Scale(cam.near_distance, cam.gaze), s);
    return ray;
}

unsigned char *RayTracer::InitializeImage(int width, int height) {
    unsigned char *image = new unsigned char[width * height * 3];
    int i = 0;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            image[i++] = 0;
            image[i++] = 0;
            image[i++] = 255;
        }
    }
    return image;
}

unsigned char *RayTracer::RenderScene(parser::Scene scene, int camera_no) {
    int width = scene.cameras.at(camera_no).image_width;
    int height = scene.cameras.at(camera_no).image_height;
    // Step 1: Initialize the image.
    unsigned char *raw_image = InitializeImage(width, height);
    // Step 2: Generate rays.
    return raw_image;
}

int main(int argc, char *argv[]) {
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
