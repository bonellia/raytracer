#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "raytracer.h"
#include "util.h"

typedef unsigned char RGB[3];

Vec3f RayTracer::Cross(const Vec3f &lhs, const Vec3f &rhs) {
    return {
            lhs.y * rhs.z - rhs.y * lhs.z,
            rhs.x * lhs.z - lhs.x * rhs.z,
            lhs.x * rhs.y - rhs.x * lhs.y};
}

float RayTracer::Dot(const Vec3f &lhs, const Vec3f &rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Vec3f RayTracer::Scale(const float &lhs, const Vec3f &rhs) {
    Vec3f vector;
    vector.x = lhs * rhs.x;
    vector.y = lhs * rhs.y;
    vector.z = lhs * rhs.z;
    return vector;
}

Vec3f RayTracer::Add(const Vec3f &lhs, const Vec3f &rhs) {
    return {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z};
}

Vec3f RayTracer::Subtract(const Vec3f &lhs, const Vec3f &rhs) {
    return {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z};
}

Vec3f RayTracer::Negate(const Vec3f &rhs) {
    return {
            -rhs.x,
            -rhs.y,
            -rhs.z};
}


float RayTracer::Length(Vec3f vector) {
    return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

Vec3f RayTracer::Normalize(Vec3f vector) {
    Vec3f normalized_vector;
    float vector_length = Length(vector);
    normalized_vector.x = vector.x / vector_length;
    normalized_vector.y = vector.y / vector_length;
    normalized_vector.z = vector.z / vector_length;
    return normalized_vector;
}

float RayTracer::Determinant(Vec3f a, Vec3f b, Vec3f c) {
    float result = 0;
    result += a.x * (a.y * b.z - b.y * a.z);
    result += a.y * (b.x * a.z - a.x * b.z);
    result += a.z * (a.x * b.y - a.y * b.x);
    return result;
}

Ray RayTracer::GenerateEyeRay(int x, int y, parser::Camera cam) {
    Ray ray;
    Vec3f su, sv, s;
    ray.origin = cam.position;
    float pixel_width = cam.near_plane.y - cam.near_plane.x / (float) cam.image_width;
    float pixel_height = (cam.near_plane.w - cam.near_plane.z) / cam.image_height;
    Vec3f cam_u = Cross(cam.gaze, cam.up);
    su = Scale(cam.near_plane.x + (x + 0.5) * (pixel_width), cam_u);
    sv = Scale(cam.near_plane.z + (x + 0.5) * (pixel_height), cam.up);
    s = Add(su, sv);

    ray.direction = Add(Scale(cam.near_distance, cam.gaze), s);
    return ray;
}

unsigned char *RayTracer::InitializeImage(int width, int height) {
    auto *image = new unsigned char[width * height * 3];
    int i = 0;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            image[i++] = 0;
            image[i++] = 255;
            image[i++] = 0;
        }
    }
    return image;
}

Touch RayTracer::IntersectSphere(Ray ray, Sphere sphere) {
    Touch touch;
    float A, B, C;
    float delta;
    Vec3f sphere_center;
    float sphere_radius;
    float t, t1, t2;

    sphere_center = this->scene.vertex_data.at(sphere.center_vertex_id);
    sphere_radius = sphere.radius;

    C += powf((ray.origin.x - sphere_center.x), 2);
    C += powf((ray.origin.y - sphere_center.y), 2);
    C += powf((ray.origin.z - sphere_center.z), 2);
    C -= powf(sphere_radius, 2);

    B += 2 * ray.direction.x * (ray.origin.x - sphere_center.x);
    B += 2 * ray.direction.y * (ray.origin.y - sphere_center.y);
    B += 2 * ray.direction.z * (ray.origin.z - sphere_center.z);

    A = powf(ray.direction.x, 2) + powf(ray.direction.y, 2) + powf(ray.direction.z, 2);

    delta = powf(B, 2) - 4 * A * C;

    if (delta > 0) {
        t1 = (-B + sqrtf(delta)) / (2 * A);
        t2 = (-B - sqrtf(delta)) / (2 * A);
        // Pick closer t value.
        t = std::min(t1, t2);
        // But this t value needs to be further than near plane.
        if (t >= 1.0) {
            // Successful hit, set touch parameters.
            touch.t = t;
            touch.position = Add(ray.origin, Scale(t, ray.direction));
            touch.normal = Normalize(Subtract(touch.position, sphere_center));
            touch.material_id = sphere.material_id;
            return touch;
        }
    }
    // If we reach here, bad luck:
    return MISS;
}

Touch RayTracer::IntersectTriangle(Ray ray, Vec3f &a, Vec3f &b, Vec3f &c, int material_id, int object_no) {
    Touch touch;
    Vec3f a_b = Subtract(a, b);
    Vec3f a_c = Subtract(a, c);
    Vec3f a_o = Subtract(b, ray.origin);

    float det_A = Determinant(a_b, a_c, ray.direction);
    if (det_A == 0.0) {
        return MISS;
    }
    float t = (Determinant(a_b, a_c, a_o)) / det_A;
    if (det_A <= 0.0) {
        return MISS;
    }
    float gamma = Determinant(a_b, a_o, ray.direction) / det_A;
    if (gamma < 0 || gamma > 1) {
        return MISS;
    }
    float beta = Determinant(a_o, a_c, ray.direction) / det_A;
    if (beta < 0 || beta > (1 - gamma)) {
        return MISS;
    }
    touch.material_id = material_id;
    touch.touched_object_no = object_no;
    touch.material_id = material_id;
    touch.t = t;
    touch.position = Add(ray.origin, Scale(t, ray.direction));
    touch.normal = Normalize(Cross(Subtract(b, a), Subtract(c, a)));
    return touch;
}

Touch RayTracer::MeshIntersection(Ray ray, int mesh_id, int material_id, int object_no) {
    Touch final_touch = MISS;
    int face_count = scene.meshes.at(mesh_id).faces.size();
    for (int face_no = 0; face_no < face_count; face_no++) {
        Touch triangle_touch;
        Vec3f v0 = scene.vertex_data[scene.meshes.at(mesh_id).faces[face_no].v0_id - 1];
        Vec3f v1 = scene.vertex_data[scene.meshes.at(mesh_id).faces[face_no].v1_id - 1];
        Vec3f v2 = scene.vertex_data[scene.meshes.at(mesh_id).faces[face_no].v2_id - 1];
        triangle_touch = IntersectTriangle(ray, v0, v1, v2, scene.meshes.at(mesh_id).material_id, object_no);
        if (triangle_touch.t >= 0) {
            if (final_touch.t == -1) {
                final_touch = triangle_touch;
            } else {
                if (final_touch.t > triangle_touch.t) {
                    final_touch = triangle_touch;
                }
            }
        }
    }
    return final_touch;
}

Vec3f RayTracer::CalculatePixelColor(Ray ray) {
    return {12, 144, 196};
}

unsigned char *RayTracer::RenderScene(Camera camera, int width, int height) {

    unsigned char *image = InitializeImage(width, height);
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            Vec3f pixel_color = {0, 0, 0};
            Ray ray = GenerateEyeRay(i, j, camera);
            pixel_color = CalculatePixelColor(ray);
            image[3 * i * width + j] = pixel_color.x;
            image[3 * i * width + j + 1] = pixel_color.y;
            image[3 * i * width + j + 2] = pixel_color.z;
        }
    }
    return image;
}

int main(int argc, char *argv[]) {
    // Sample usage for reading an XML scene file
    Scene scene;
    util::Util util;
    RayTracer ray_tracer;
    scene.loadFromXml(argv[1]);
    ray_tracer.scene = scene;
    util.PrintSceneDetails(scene);
    // TODO: Add multithreading.
    for (auto camera: scene.cameras) {
        auto width = camera.image_width;
        auto height = camera.image_height;
        unsigned char *image = ray_tracer.RenderScene(camera, width, height);
        write_ppm(camera.image_name.c_str(), image, width, height);
    }
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

}
