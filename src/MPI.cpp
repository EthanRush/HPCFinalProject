//==============================================================================================
// Based on Code Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "rtweekend.h"

#include "box.h"
#include "bvh.h"
#include "camera.h"
#include "color.h"
#include "constant_medium.h"
#include "hittable_list.h"
#include "material.h"
#include "moving_sphere.h"
#include "sphere.h"
#include "texture.h"

#include <iostream>
#include <fstream>

#include <mpi.h>
#include <time.h>
#include <string>
#include <cstring>

double CLOCK() {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (t.tv_sec * 1000) + (t.tv_nsec * 1e-6);
}


color ray_color(const ray& r, const color& background, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0,0,0);

    // If the ray hits nothing, return the background color.
    if (!world.hit(r, 0.001, infinity, rec))
        return background;

    ray scattered;
    color attenuation;
    color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    return emitted + attenuation * ray_color(scattered, background, world, depth-1);
}


hittable_list random_scene() {
    hittable_list world;

    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));

    world.add(make_shared<sphere>(point3(0,-1000,0), 1000, make_shared<lambertian>(checker)));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

            if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    auto center2 = center + vec3(0, random_double(0,.5), 0);
                    world.add(make_shared<moving_sphere>(
                        center, center2, 0.0, 1.0, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return hittable_list(make_shared<bvh_node>(world, 0.0, 1.0));
}


hittable_list cornell_box() {
    hittable_list objects;

    auto red   = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light = make_shared<diffuse_light>(color(15, 15, 15));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<xz_rect>(213, 343, 227, 332, 554, light));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    shared_ptr<hittable> box1 = make_shared<box>(point3(0,0,0), point3(165,330,165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265,0,295));
    objects.add(box1);

    shared_ptr<hittable> box2 = make_shared<box>(point3(0,0,0), point3(165,165,165), white);
    box2 = make_shared<rotate_y>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(130,0,65));
    objects.add(box2);

    return objects;
}


hittable_list final_scene() {
    hittable_list boxes1;
    auto ground = make_shared<lambertian>(color(0.48, 0.83, 0.53));

    const int boxes_per_side = 20;
    for (int i = 0; i < boxes_per_side; i++) {
        for (int j = 0; j < boxes_per_side; j++) {
            auto w = 100.0;
            auto x0 = -1000.0 + i*w;
            auto z0 = -1000.0 + j*w;
            auto y0 = 0.0;
            auto x1 = x0 + w;
            auto y1 = random_double(1,101);
            auto z1 = z0 + w;

            boxes1.add(make_shared<box>(point3(x0,y0,z0), point3(x1,y1,z1), ground));
        }
    }

    hittable_list objects;

    objects.add(make_shared<bvh_node>(boxes1, 0, 1));

    auto light = make_shared<diffuse_light>(color(7, 7, 7));
    objects.add(make_shared<xz_rect>(123, 423, 147, 412, 554, light));

    auto center1 = point3(400, 400, 200);
    auto center2 = center1 + vec3(30,0,0);
    auto moving_sphere_material = make_shared<lambertian>(color(0.7, 0.3, 0.1));
    objects.add(make_shared<moving_sphere>(center1, center2, 0, 1, 50, moving_sphere_material));

    objects.add(make_shared<sphere>(point3(260, 150, 45), 50, make_shared<dielectric>(1.5)));
    objects.add(make_shared<sphere>(
        point3(0, 150, 145), 50, make_shared<metal>(color(0.8, 0.8, 0.9), 1.0)
    ));

    auto boundary = make_shared<sphere>(point3(360,150,145), 70, make_shared<dielectric>(1.5));
    objects.add(boundary);
    objects.add(make_shared<constant_medium>(boundary, 0.2, color(0.2, 0.4, 0.9)));
    boundary = make_shared<sphere>(point3(0,0,0), 5000, make_shared<dielectric>(1.5));
    objects.add(make_shared<constant_medium>(boundary, .0001, color(1,1,1)));

    auto pertext = make_shared<noise_texture>(0.1);
    objects.add(make_shared<sphere>(point3(220,280,300), 80, make_shared<lambertian>(pertext)));

    hittable_list boxes2;
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    int ns = 1000;
    for (int j = 0; j < ns; j++) {
        boxes2.add(make_shared<sphere>(point3::random(0,165), 10, white));
    }

    objects.add(make_shared<translate>(
        make_shared<rotate_y>(
            make_shared<bvh_node>(boxes2, 0.0, 1.0), 15),
            vec3(-100,270,395)
        )
    );

    return objects;
}


void render(std::ostream &out, hittable_list world, camera cam, float aspect_ratio, int image_width, int samples_per_pixel , int max_depth , color background) {
   
    int world_rank;
    int world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size != 12) {
        fprintf(stderr, "World size must be twelve");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    
    double start, finish, total;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    if (world_size != 12) {
        fprintf(stderr, "World size must be twelve");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (world_rank == 0) {
        start = CLOCK();
        out << "P3\n" << image_width << ' ' << image_height << "\n255\n";
    }

    int rows_per_proc = int(std::ceil(image_height * 1.0 / world_size));
    int start_index, end_index;


    // If image not divisible by world size need to allocate the extra rows
    int remainder = image_height % world_size;
    if (remainder != 0) {
         start_index = (image_height - 1) - ((rows_per_proc * world_rank) - std::max(0, (world_rank - remainder)));
        // height 50: 0 = 49, 1 = 44, 2 = 39, 3 = 35, 4 = 31, 5 = 27, 6 = 23, 7 = 19, 8 = 15, 9 = 11, 10 = 7, 11 = 3
        
        end_index = (image_height - 1) - (((rows_per_proc * (world_rank + 1)) - std::max(0, ((world_rank + 1) - remainder))) - 1) ;
        // height 50:  0 = 45, 1 = 40, 2 = 36, 3 = 32, ..., 11 = 0
    }
    else {
        start_index = (image_height - 1) - (rows_per_proc * world_rank);
        // height 24: 0 = 23, 1 = 21, 2 = 19, 3 = 17, ..., 11 = 1
        end_index = (image_height - 1) - (rows_per_proc * (world_rank + 1) - 1);
        // height 24: 0 = 22, 1 = 20, 2 = 18, ..., 11 = 0
    }
    std::string localstr = "";

    for (int j = start_index; j >= end_index; --j) {
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, background, world, max_depth);
            }
            write_color(localstr, pixel_color, samples_per_pixel);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int mylen = localstr.length();
  
    int* recvcounts = NULL;

    /* Only root has the received data */
    if (world_rank == 0) {
        recvcounts =(int *) malloc(world_size * sizeof(int));
    }

    // Gathers length of each process' string for displacement calculations
    MPI_Gather(&mylen, 1, MPI_INT,
        recvcounts, 1, MPI_INT,
        0, MPI_COMM_WORLD);
    
    int totlen = 0;
    int* displs = NULL;
    char* totalstring = NULL;

    if (world_rank == 0) {
        displs = (int *)malloc(world_size * sizeof(int));

        displs[0] = 0;
        totlen += recvcounts[0] + 1;

        for (int i = 1; i < world_size; i++) {
            totlen += recvcounts[i] + 1;   /* plus one for space or \0 after words */
            displs[i] = displs[i - 1] + recvcounts[i - 1] + 1;
        }

        /* allocate string, pre-fill with spaces and null terminator */
        totalstring = (char *)malloc(totlen * sizeof(char));
        for (int i = 0; i < totlen - 1; i++)
            totalstring[i] = ' ';
        totalstring[totlen - 1] = '\0';
    }
    // Convert the string to a char * for use with Gather
    const char* localCharArr = localstr.c_str();

    MPI_Gatherv(localCharArr, mylen, MPI_CHAR,
        totalstring, recvcounts, displs, MPI_CHAR,
        0, MPI_COMM_WORLD);


    MPI_Finalize();
    if (world_rank == 0) {
        finish = CLOCK();
        total = finish - start;
        std::cout << "Total Render Time: " << total << std::endl;
        std::cout << "Begin write to file" << std::endl;
        out << totalstring;
        std::cout << "Finished write to file" << std::endl;
    }

}

int main(int argc, char* argv[]) {

    if (argc > 2) {
        std::cout << "Too many arguments, please just enter the scene number (1-3) or nothing to default to 1" << std::endl;
        return -1;
    }

    int sceneNum = std::stoi(argv[1]);

    int world_rank;
    int world_size;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    unsigned int seed;
    // Ensure that all processes are using same seed for their randomness
    if (world_rank == 0) {
        seed = time(NULL);
    }
    MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    srand(seed);
    // Image

    auto aspect_ratio = 16.0 / 9.0;
    int image_width = 400;
    int samples_per_pixel = 100;
    int max_depth = 50;

    // World

    hittable_list world;

    // Camera
    const vec3 vup(0, 1, 0);
    const auto dist_to_focus = 10.0;

    point3 lookfrom;
    point3 lookat;
    auto vfov = 40.0;
    auto aperture = 0.0;
    color background(0,0,0);
    std::ofstream out;
    std::string fileName;
   
    switch (sceneNum)
    {
        case 1:
        default:
        {
            // Scene 1
            world = random_scene();
            background = color(0.70, 0.80, 1.00);
            lookfrom = point3(13, 2, 3);
            lookat = point3(0, 0, 0);
            vfov = 20.0;
            aperture = 0.1;
            fileName = "shirleyOrbs_mpi.ppm";
            break;
        }
        case 2:
        {
            // Scene 2
            world = cornell_box();
            aspect_ratio = 1.0;
            image_width = 600;
            samples_per_pixel = 200;
            lookfrom = point3(278, 278, -800);
            lookat = point3(278, 278, 0);
            vfov = 40.0;
            fileName = "cornell_mpi.ppm";
            break;
        }
        case 3:
        {
            // Scene 3
            world = final_scene();
            aspect_ratio = 1.0;
            image_width = 800;
            samples_per_pixel = 600;
            lookfrom = point3(478, 278, -600);
            lookat = point3(278, 278, 0);
            vfov = 40.0;
            fileName = "final_mpi.ppm";
            break;
        }
    }
    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

    if (world_rank == 0) {
        out.open(fileName);
    }

    render(out, world, cam, aspect_ratio, image_width, samples_per_pixel, max_depth, background);
    
    if (world_rank == 0) {
        out.close();
    }
   
    return 0;
}
