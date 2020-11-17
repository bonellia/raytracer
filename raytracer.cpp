#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "raytracer.h"
#include "util.h"
#include <thread>
#include <future>

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
    result.y = lhs.y * rhs.y;
    result.z = lhs.z * rhs.z;
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
    normalized_vector = Circumsize(vector, vector_length);
    return normalized_vector;
}

float RayTracer::Determinant(const Vec3f &a, const Vec3f &b, const Vec3f &c) {
    float result = a.x * (b.y * c.z - c.y * b.z) + a.y * (c.x * b.z - b.x * c.z) + a.z * (b.x * c.y - b.y * c.x);
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
            image[i++] = 0;
            image[i++] = 0;
        }
    }
    return image;
}

TouchAttempt RayTracer::SphereIntersectionTest(const Ray &ray, const Sphere &sphere) {
    TouchAttempt touch_attempt;
    Vec3f sphere_center = this->scene.vertex_data.at(sphere.center_vertex_id - 1);
    const float sphere_radius = sphere.radius;
    const Vec3f origin_minus_center = Subtract(ray.origin, sphere_center);
    const float A = Dot(ray.direction, ray.direction);
    const float B = 2 * Dot(ray.direction, origin_minus_center);
    const float C = Dot(origin_minus_center, origin_minus_center) - sphere_radius * sphere_radius;
    const float delta = B * B - 4 * A * C;

    if (delta > 0) {
        // Successful hit, set touch attempt parameters.
        // Pick closer t value.
        const float t1 = (-1 * B + std::sqrt(delta)) / 2 * A;
        const float t2 = (-1 * B - std::sqrt(delta)) / 2 * A;
        float t = std::min(t1, t2);
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

TouchAttempt RayTracer::TriangleIntersectionTest(const Ray &ray, const Vec3f &a, const Vec3f &b, const Vec3f &c,
                                                 const int material_id) {
    TouchAttempt touch_attempt = MISS;
    Vec3f a_b = Subtract(a, b);
    Vec3f a_c = Subtract(a, c);
    Vec3f a_o = Subtract(a, ray.origin);

    float det_A = Determinant(a_b, a_c, ray.direction);
    if (det_A == 0.0) {
        return MISS;
    }
    float t = (Determinant(a_b, a_c, a_o)) / det_A;
    if (t <= 0.0) {
        return MISS;
    }
    float gamma = Determinant(a_b, a_o, ray.direction) / det_A;
    if (gamma < 0 || gamma > 1) {
        return MISS;
    }
    float beta = Determinant(a_o, a_c, ray.direction) / det_A;
    if (beta < 0 || beta > 1 || beta + gamma > 1) {
        return MISS;
    }
    touch_attempt.material_id = material_id;
    touch_attempt.t = t;
    touch_attempt.position = Add(ray.origin, Scale(t, ray.direction));
    touch_attempt.normal = Normalize(Cross(Subtract(b, a), Subtract(c, a)));
    touch_attempt.contact = BOOB;
    return touch_attempt;
}

TouchAttempt RayTracer::MeshIntersectionTest(const Ray &ray, const Mesh &mesh) {
    TouchAttempt closest_touch_attempt = MISS;
    closest_touch_attempt.t = std::numeric_limits<float>::max();
    int face_count = mesh.faces.size();
    for (int face_no = 0; face_no < face_count; face_no++) {
        TouchAttempt current_attempt;
        Vec3f v0 = scene.vertex_data[mesh.faces[face_no].v0_id - 1];
        Vec3f v1 = scene.vertex_data[mesh.faces[face_no].v1_id - 1];
        Vec3f v2 = scene.vertex_data[mesh.faces[face_no].v2_id - 1];
        current_attempt = TriangleIntersectionTest(ray, v0, v1, v2, mesh.material_id);
        if (current_attempt.contact == BOOB && current_attempt.t < closest_touch_attempt.t) {
            closest_touch_attempt = current_attempt;
            closest_touch_attempt.contact = BODY;
        } else {
            continue;
        }
    }
    return closest_touch_attempt;
}

TouchAttempt RayTracer::FindClosestContact(const Ray &ray) {
    TouchAttempt closest_touch_attempt = MISS;
    closest_touch_attempt.t = std::numeric_limits<float>::max();
    // Check spheres.
    for (const auto &sphere : scene.spheres) {
        TouchAttempt current_attempt = SphereIntersectionTest(ray, sphere);
        if (current_attempt.contact == BUM && current_attempt.t < closest_touch_attempt.t) {
            closest_touch_attempt = current_attempt;
        } else {
            continue;
        }
    }

    // Check triangles.
    for (const auto &triangle : scene.triangles) {
        Vec3f a = scene.vertex_data[triangle.indices.v0_id - 1];
        Vec3f b = scene.vertex_data[triangle.indices.v1_id - 1];
        Vec3f c = scene.vertex_data[triangle.indices.v2_id - 1];
        TouchAttempt current_attempt = TriangleIntersectionTest(ray, a, b, c, triangle.material_id);
        if (current_attempt.contact == BOOB && current_attempt.t < closest_touch_attempt.t) {
            closest_touch_attempt = current_attempt;
        } else {
            continue;
        }
    }
    // Check meshes.
    for (const auto &mesh : scene.meshes) {
        TouchAttempt current_attempt = MeshIntersectionTest(ray, mesh);
        if (current_attempt.contact == BODY && current_attempt.t < closest_touch_attempt.t) {
            closest_touch_attempt = current_attempt;
        } else {
            continue;
        }
    }
    return closest_touch_attempt;
}

Vec3f RayTracer::CalculatePixelColor(const Ray &ray, const TouchAttempt &touch_attempt, const Camera &cam, int depth) {
    Vec3f pixel_color = {0, 0, 0};
    if (depth > scene.max_recursion_depth)
        return pixel_color;
    if (touch_attempt.contact != NAH) {
        Material touched_mat = scene.materials[touch_attempt.material_id - 1];
        // Handle reflections first.
        if (Length(touched_mat.mirror) > 0) {
            Vec3f reflection_vector = Subtract(ray.direction,
                                               Scale(
                                                       2.0f * (Dot(touch_attempt.normal, ray.direction)),
                                                       touch_attempt.normal));
            Ray reflection_ray;
            reflection_ray.origin = Add(touch_attempt.position, Scale(scene.shadow_ray_epsilon, reflection_vector));
            reflection_ray.direction = Normalize(reflection_vector);
            Vec3f mirror_component = VectorScale(touched_mat.mirror,
                                                 CalculatePixelColor(reflection_ray, touch_attempt, cam, depth + 1));

            pixel_color = Add(pixel_color, mirror_component);
        }

        // Ambient shading next.
        Vec3f ambient_component = VectorScale(scene.ambient_light, touched_mat.ambient);
        pixel_color = Add(pixel_color, ambient_component);
        // Then lighting.
        for (auto light : scene.point_lights) {
            Vec3f light_vector = Subtract(light.position, touch_attempt.position);
            float light_distance = Length(light_vector);
            Vec3f irradiance = Circumsize(light.intensity,
                                          powf(Length(Subtract(touch_attempt.position, light.position)), 2));
            Ray light_ray;
            light_ray.origin = Add(touch_attempt.position, Scale(scene.shadow_ray_epsilon, Normalize(light_vector)));
            light_ray.direction = Normalize(light_vector);

            // Shadows.
            TouchAttempt shadow_touch = FindClosestContact(light_ray);
            if (shadow_touch.t > 0 && shadow_touch.t <= light_distance) {
                continue;
            }

            // Diffuse.
            Vec3f diffuse_component = Scale(std::max(0.0f, Dot(touch_attempt.normal, Normalize(light_vector))),
                                            VectorScale(touched_mat.diffuse, irradiance));
            pixel_color = Add(pixel_color, diffuse_component);

            // Specular (Blinn-Phong).
            Vec3f h = Add(Negate(Normalize(ray.direction)), Normalize(light_vector));
            h = Normalize(h);
            Vec3f specular_component = Scale(
                    powf(std::max(0.0f, Dot(touch_attempt.normal, h)), touched_mat.phong_exponent),
                    VectorScale(touched_mat.specular, irradiance));
            pixel_color = Add(pixel_color, specular_component);
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

void RayTracer::SetImagePixelRGB(unsigned char *&image, int row, int column, const int width, const RGB color) {
    image[3 * (row * width + column)] = color[0];
    image[3 * (row * width + column) + 1] = color[1];
    image[3 * (row * width + column) + 2] = color[2];
}

void RayTracer::RenderBar(const Camera &cam, unsigned char *&image, const int width_from, const int width_to,
                          const int height) {
    for (int row = 0; row < height; row++) {
        for (int column = width_from; column < width_to; column++) {
            Vec3f pixel_color;
            Ray ray = GenerateEyeRay(row, column, cam);
            TouchAttempt touch_attempt = MISS;
            touch_attempt = FindClosestContact(ray);
            pixel_color = CalculatePixelColor(ray, touch_attempt, cam, 0);
            RGB clamped_pixel_color;
            clamped_pixel_color[0] = pixel_color.x > 255 ? 255 : (int) ceilf(pixel_color.x);
            clamped_pixel_color[1] = pixel_color.y > 255 ? 255 : (int) ceilf(pixel_color.y);
            clamped_pixel_color[2] = pixel_color.z > 255 ? 255 : (int) ceilf(pixel_color.z);

            /*
            if (clamped_pixel_color[0] != 0 || clamped_pixel_color[1] != 0 || clamped_pixel_color[2] != 0) {
                std::cout << clamped_pixel_color[0]<< ' ' << clamped_pixel_color[1] << ' ' << clamped_pixel_color[2] << ' ' << std::endl;
                std::cout << row << ' ' << column << std::endl;
            }
            */
            SetImagePixelRGB(image, row, column, cam.image_width, clamped_pixel_color);
        }
    }
}

unsigned char *RayTracer::RenderScene(const Camera &cam, const int width, const int height) {

    unsigned char *image = InitializeImage(width, height);

    const unsigned int thread_count = std::max(std::thread::hardware_concurrency(), 1u);
    std::cout << "Thread count is: " << thread_count << std::endl;
    const int stride = static_cast<int>(width / thread_count);
    {
        std::vector<std::future<void>> partial_render_tasks;

        for (std::size_t i = 0; i < thread_count; ++i) {
            partial_render_tasks.push_back(std::async(std::launch::async, [&, i] {
                RenderBar(cam, image, i * stride, (i + 1) * stride, height);
            }));
        }
        // For a potential leftover bar in case image width is not a multitude of thread_count.
        if (width % thread_count) {
            partial_render_tasks.push_back(std::async(std::launch::async, [&] {
                RenderBar(cam, image, thread_count * stride, width, height);
            }));
        }
    }
    return image;
}

int main(int argc, char *argv[]) {
    // Sample usage for reading an XML scene file
    // Better make use of argc:
    if (argc < 2) {
        printf("Usage: ./raytracer <scene file>\n");
        return 1;
    }
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
