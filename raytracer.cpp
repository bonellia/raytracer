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
    Vec3f result;
    result.x = lhs * rhs.x;
    result.y = lhs * rhs.y;
    result.z = lhs * rhs.z;
    return result;
}

Vec3f RayTracer::VectorScale(const Vec3f &lhs, const Vec3f &rhs) {
    Vec3f result;
    result.x = lhs.x * rhs.x;
    result.y = lhs.x * rhs.y;
    result.z = lhs.x * rhs.z;
    return result;
}

Vec3f RayTracer::Circumsize(const Vec3f &lhs, const float &rhs) {
    Vec3f result;
    result.x = lhs.x / rhs;
    result.y = lhs.y / rhs;
    result.z = lhs.z / rhs;
    return result;
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

Vec3f RayTracer::Negate(const Vec3f &vector) {
    Vec3f result;
    result.x = -vector.x;
    result.y = -vector.y;
    result.z = -vector.z;
    return result;
}


float RayTracer::Length(const Vec3f &vector) {
    return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

Vec3f RayTracer::Normalize(const Vec3f &vector) {
    Vec3f normalized_vector;
    float vector_length = Length(vector);
    normalized_vector.x = vector.x / vector_length;
    normalized_vector.y = vector.y / vector_length;
    normalized_vector.z = vector.z / vector_length;
    return normalized_vector;
}

float RayTracer::Determinant(const Vec3f &a, const Vec3f &b, const Vec3f &c) {
    float result = 0;
    result += a.x * (a.y * b.z - b.y * a.z);
    result += a.y * (b.x * a.z - a.x * b.z);
    result += a.z * (a.x * b.y - a.y * b.x);
    return result;
}

Ray RayTracer::GenerateEyeRay(int pixel_row, int pixel_column, const parser::Camera &cam) {
    // We find the position of the direction vector by traversing near plane.
    Ray ray;
    Vec3f plane_center = Add(cam.position, Scale(cam.near_distance, cam.gaze));
    Vec3f cam_right_vector = Normalize(Cross(cam.gaze, cam.up));
    Vec3f plane_top_left = Add(plane_center,
                               Add(Scale(cam.near_plane.x, cam_right_vector), Scale(cam.near_plane.w, cam.up)));
    float u = (float) (pixel_column + 0.5) * (cam.near_plane.y - cam.near_plane.x) / (float) cam.image_width;
    float v = (float) (pixel_row + 0.5) * (cam.near_plane.w - cam.near_plane.z) / (float) cam.image_height;

    Vec3f pixel_position = Add(plane_top_left, Subtract(Scale(u, cam_right_vector), Scale(v, cam.up)));
    ray.origin = cam.position;
    ray.direction = Normalize(Subtract(pixel_position, ray.origin));
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

TouchAttempt RayTracer::SphereIntersectionTest(const Ray &ray, const Sphere &sphere) {
    TouchAttempt touch_attempt;
    float delta = 0;
    Vec3f sphere_center = this->scene.vertex_data.at(sphere.center_vertex_id - 1);
    const float sphere_radius = sphere.radius;
    const Vec3f origin_minus_center = Subtract(ray.origin, sphere_center);
    float t;
    const float A = Dot(ray.direction, ray.direction);
    const float B = 2 * Dot(ray.direction, origin_minus_center);
    const float C = Dot(origin_minus_center, origin_minus_center) - sphere_radius * sphere_radius;
    delta = B * B - 4 * A * C;

    if (delta > 0) {
        // Successful hit, set touch attempt parameters.
        // Pick closer t value.
        const float t1 = (-1 * B + std::sqrt(delta)) / 2 * A;
        const float t2 = (-1 * B - std::sqrt(delta)) / 2 * A;
        t = std::min(t1, t2);
        touch_attempt.t = t;
        touch_attempt.position = Add(ray.origin, Scale(t, ray.direction));
        touch_attempt.normal = Circumsize(Subtract(touch_attempt.position, sphere_center), sphere_radius);
        touch_attempt.material_id = sphere.material_id;
        touch_attempt.contact = BUM;
        return touch_attempt;
    }
    // If we reach here, bad luck:
    return MISS;
}

TouchAttempt RayTracer::TriangleIntersectionTest(const Ray ray, const Vec3f &a, const Vec3f &b, const Vec3f &c,
                                                 const int material_id) {
    TouchAttempt touch;
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
    touch.material_id = material_id;
    touch.t = t;
    touch.position = Add(ray.origin, Scale(t, ray.direction));
    touch.normal = Normalize(Cross(Subtract(b, a), Subtract(c, a)));
    return touch;
}

TouchAttempt RayTracer::MeshIntersectionTest(const Ray &ray, const Mesh &mesh) {
    TouchAttempt touch_attempt = MISS;
    int face_count = mesh.faces.size();
    for (int face_no = 0; face_no < face_count; face_no++) {
        TouchAttempt triangle_touch;
        Vec3f v0 = scene.vertex_data[mesh.faces[face_no].v0_id - 1];
        Vec3f v1 = scene.vertex_data[mesh.faces[face_no].v1_id - 1];
        Vec3f v2 = scene.vertex_data[mesh.faces[face_no].v2_id - 1];
        triangle_touch = TriangleIntersectionTest(ray, v0, v1, v2, mesh.material_id);
        if (triangle_touch.t >= 0) {
            if (touch_attempt.t == -1) {
                continue;
            } else {
                if (touch_attempt.t > triangle_touch.t) {
                    touch_attempt = triangle_touch;
                }
            }
        }
    }
    return touch_attempt;
}

TouchAttempt RayTracer::FindClosestContact(const Ray &ray) {
    TouchAttempt closest_touch_attempt = MISS;
    // Check spheres.
    for (const auto &sphere : scene.spheres) {
        TouchAttempt current_attempt = SphereIntersectionTest(ray, sphere);
        if (current_attempt.contact == BUM && current_attempt.t < closest_touch_attempt.t) {
            closest_touch_attempt = current_attempt;
        } else {
            continue;
        }
    }
    // Check meshes.
    for (const auto &mesh : scene.meshes) {
        TouchAttempt current_attempt = MeshIntersectionTest(ray, mesh);
        if (current_attempt.contact == BUM && current_attempt.t < closest_touch_attempt.t) {
            closest_touch_attempt = current_attempt;
        } else {
            continue;
        }
    }
    // Check triangles.
    for (const auto &triangle : scene.triangles) {
        Vec3f v0 = scene.vertex_data[triangle.indices.v0_id - 1];
        Vec3f v1 = scene.vertex_data[triangle.indices.v1_id - 1];
        Vec3f v2 = scene.vertex_data[triangle.indices.v2_id - 1];
        TouchAttempt current_attempt = TriangleIntersectionTest(ray, v0, v1, v2, triangle.material_id);
        if (current_attempt.contact == BUM && current_attempt.t < closest_touch_attempt.t) {
            closest_touch_attempt = current_attempt;
        } else {
            continue;
        }
    }
    return closest_touch_attempt;
}

Vec3f RayTracer::CalculatePixelColor(const Ray &ray, int depth) {
    Vec3f pixel_color = {0, 0, 0};
    TouchAttempt closest_touch = MISS;
    closest_touch = FindClosestContact(ray);
    if (closest_touch.t < 0) {
        Material touched_mat = scene.materials[closest_touch.material_id - 1];

        // Handle reflections early.
        if (Length(touched_mat.mirror) > 0) {
            Vec3f reflection_vector = Subtract(ray.direction,
                                               Scale(
                                                       2.0f * (Dot(closest_touch.normal, ray.direction)),
                                                       closest_touch.normal));
            Ray reflection_ray;
            reflection_ray.origin = Add(closest_touch.position, Scale(scene.shadow_ray_epsilon, reflection_vector));
            reflection_ray.direction = Normalize(reflection_vector);
            Vec3f mirror_component = VectorScale(touched_mat.mirror, CalculatePixelColor(reflection_ray, depth + 1));

            pixel_color = Add(pixel_color, mirror_component);
        }
        // Ambient shading after.
        Vec3f ambient_component = VectorScale(scene.ambient_light, touched_mat.ambient);
        pixel_color = Add(pixel_color, ambient_component);
        // Then lighting.
        for (auto light : scene.point_lights) {
            Vec3f light_vector = Subtract(light.position, closest_touch.position);
            float light_distance = Length(light_vector);
            Vec3f light_component = Circumsize(light.intensity,
                                               powf(Length(Subtract(closest_touch.position, light.position)), 2));
            Ray light_ray;
            light_ray.origin = Add(closest_touch.position, Scale(scene.shadow_ray_epsilon, Normalize(light_vector)));
            light_ray.direction = Normalize(light_vector);

            // Shadows.
            bool skip_rest = false;
            TouchAttempt shadow_touch = FindClosestContact(light_ray);
            if (shadow_touch.t > 0 && shadow_touch.t <= light_distance) {
                continue;
            }

            // Diffuse.
            Vec3f diffuse_component = Scale(std::max(0.0f, Dot(closest_touch.normal, Normalize(light_vector))),
                                            VectorScale(touched_mat.diffuse, light_component));
            pixel_color = Add(pixel_color, diffuse_component);

            // Specular (Blinn-Phong).
            Vec3f h = Add(Negate(Normalize(ray.direction)), Normalize(light_vector));
            h = Circumsize(h, Length(h));
            Vec3f specular_component = Add(pixel_color, specular_component);
        }
        pixel_color.x = pixel_color.x > 255 ? 255 : pixel_color.x;
        pixel_color.y = pixel_color.y > 255 ? 255 : pixel_color.y;
        pixel_color.z = pixel_color.z > 255 ? 255 : pixel_color.z;
        return pixel_color;
    }
    pixel_color.x = (float) scene.background_color.x;
    pixel_color.y = (float) scene.background_color.y;
    pixel_color.z = (float) scene.background_color.z;
    return pixel_color;
}

unsigned char *RayTracer::RenderScene(const Camera &camera, const int width, const int height) {

    unsigned char *image = InitializeImage(width, height);
    for (int row = 0; row < height; row++) {
        for (int column = 0; column < width; column++) {
            Vec3f pixel_color;
            Ray ray = GenerateEyeRay(row, column, camera);
            pixel_color = CalculatePixelColor(ray, 0);
            image[3 * (row * width + column)] = pixel_color.x;
            image[3 * (row * width + column) + 1] = pixel_color.y;
            image[3 * (row * width + column) + 2] = pixel_color.z;
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
