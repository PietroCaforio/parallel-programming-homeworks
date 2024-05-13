#include <unistd.h>
#include <string.h>
#include "raytracer.h"
#include <atomic>

//BLAME:Pietro
#define NUM_THREADS 50
#define NUM_JOBS 600
#define NUM_SAMPLES2 5

inline void set_hit_info(float root, const Ray& ray, const Sphere& sphere, Hit& hit) {
    hit.t = root;
    hit.point = ray_at(ray, hit.t);
    Vector3 outward_normal = (hit.point - sphere.center) / sphere.radius;
    set_face_normal(hit, ray, outward_normal);
}
bool sphere_hit2(const Sphere& sphere, const Ray& ray, float t_min, float t_max, Hit& hit) {
    Vector3 diff = ray.origin_point - sphere.center;

    float a = dot(ray.direction, ray.direction);
    float b = dot(diff, ray.direction);
    float c = dot(diff, diff) - sphere.radius * sphere.radius;

    float discriminant = b * b - a * c;

    if (discriminant >= 0) {
        float sqrt_discriminant = std::sqrt(discriminant);
        float inv_a = 1.0f / a;

        float root = (-b - sqrt_discriminant) * inv_a;
        if (root > t_min && root < t_max) {
            set_hit_info(root, ray, sphere, hit);
            return true;
        }

        root = (-b + sqrt_discriminant) * inv_a;
        if (root > t_min && root < t_max) {
            set_hit_info(root, ray, sphere, hit);
            return true;
        }
    }
    return false;
}

/*
** Checks if the given ray hits a sphere surface and returns.
** Also returns hit data which contains material information.
*/
bool check_sphere_hit(const std::vector<Sphere>& spheres, const Ray& ray, float t_min, float t_max, Hit& hit) {
    Hit closest_hit;
    bool has_hit = false;
    auto closest_hit_distance = t_max;
    Material material;

    for(std::size_t i = 0; i < spheres.size(); i++) {
        const auto& sphere = spheres[i];
        if(sphere_hit2(sphere, ray, t_min, closest_hit_distance, closest_hit)) {
            has_hit = true;
            closest_hit_distance = closest_hit.t;
            material = sphere.material;
        }
    }

    if(has_hit) {
        hit = closest_hit;
        hit.material = material;
    }
    
    return has_hit;
}

/*
** Traces a ray, returns color for the corresponding pixel.
*/
Vector3 trace_ray(const Ray& ray, const std::vector<Sphere>& spheres, int depth) {
    if (depth <= 0) {
        return Vector3(0, 0, 0);
    }

    Hit hit;
    if(check_sphere_hit(spheres, ray, 0.002f, FLT_MAX, hit)) {
        Ray outgoing_ray;
        Vector3 attenuation;

        if(metal_scater(hit.material, ray, hit, attenuation, outgoing_ray)) {
            auto ray_color = trace_ray(outgoing_ray, spheres, depth - 1);
            return Vector3(ray_color.x * attenuation.x, ray_color.y * attenuation.y, ray_color.z * attenuation.z);
        }

        return Vector3(0, 0, 0);
    }

    Vector3 unit_direction = unit_vector(ray.direction);
    auto t = 0.5 * (unit_direction.y + 1.0);
    return Vector3(1.0, 1.0, 1.0) * (1.0 - t) + Vector3(0.5, 0.75, 1.0) * t;
}

//-------------------------Parallelization functions and data structures -----------------------------------

struct ThreadData {
    int thread_id;
    int height;
    int width;
    int samples;
    int depth;
    Camera* camera;
    std::vector<Sphere>* spheres;
    int* image_data; //shared image data
    Checksum* checksum; //For now we leave it like that
    //pthread_mutex_t* mutex;
    ThreadData(int thread_id, int height,int width,int samples,int depth,Camera* camera,std::vector<Sphere>* spheres,int* image_data, Checksum* checksum)
    :thread_id(thread_id),height(height),width(width),samples(samples),depth(depth),camera(camera),spheres(spheres),image_data(image_data),checksum(checksum){

    }
    ThreadData(){

    }
    
};

std::atomic <int> next_job = 0;


void* perform_work(void* argument){
    while (true) {
        int thread_id = next_job++;
        if (thread_id >= NUM_JOBS) {break;}
        ThreadData* data = (struct ThreadData*)argument;
        //int thread_id = data->thread_id;
        int height = data->height;
        int width = data->width;
        int samples = data->samples;
        int depth = data->depth;
        Camera camera = *(data -> camera);
        std::vector<Sphere>& spheres = *(data->spheres);
        int* image_data = data->image_data;
        Checksum& checksum = *(data->checksum);
        
        int interval = height/NUM_JOBS; 
        for(int y = height - 1 - thread_id*interval; y >= height - (thread_id+1)*interval; y--) {
            for(int x = 0; x < width; x++) {
                Vector3 pixel_color(0,0,0);
                for(int s = 0; s < samples; s++) {
                    auto u = (float) (x + random_float()) / (width - 1);
                    auto v = (float) (y + random_float()) / (height - 1);
                    auto r = get_camera_ray(camera, u, v);
                    pixel_color += trace_ray(r, spheres, depth);
                }
                auto output_color = compute_color(checksum, pixel_color, samples);
                int pos = ((height - 1 - y) * width + x) * 3;
                image_data[pos] = output_color.r;
                image_data[pos + 1] = output_color.g;
                image_data[pos + 2] = output_color.b;
            }
        }
    }
    
    pthread_exit(NULL);
}

int main(int argc, char **argv) {
    int width = IMAGE_WIDTH;
    int height = IMAGE_HEIGHT;
    int samples = NUM_SAMPLES2;
    int depth = SAMPLE_DEPTH;

    // This option parsing is not very interesting.
    int no_output = 0;
    char file_name[256] = "render.ppm";
    int c;
    while ((c = getopt(argc, argv, "d:s:r:n:f:")) != -1)
    {
        switch (c)
        {
            case 'd':
                if (sscanf(optarg, "%d", &depth) != 1)
                    goto error;
                break;
            case 's':
                if (sscanf(optarg, "%d", &samples) != 1)
                    goto error;
                break;
            case 'r':
                if (sscanf(optarg, "%dx%d", &width, &height) != 2)
                    goto error;
                break;
            case 'n':
                if (sscanf(optarg, "%d", &no_output) != 1)
                    goto error;
                break;
            case 'f':
                strncpy(file_name, optarg, sizeof(file_name));
                file_name[255] = '\0'; // safe-guard null-terminator to disable gcc warning
                break;
            case '?':
            error: fprintf(stderr,
                        "Usage:\n"
                        "-d \t number of times a ray can bounce\n"
                        "-s \t number of samples per pixel\n"
                        "-r \t image resolution to be computed\n"
                        "-f \t output file name\n"
                        "-n \t no output(default: 0)\n"
                        "\n"
                        "Example:\n"
                        "%s -d 10 -s 50 -r 720x480 -f tracer.ppm\n",
                        argv[0]);
                exit(EXIT_FAILURE);
                break;
        }
    }

    // Calculating the aspect ratio and creating the camera for the rendering
    const auto aspect_ratio = (float) width / height;
    Camera camera(Vector3(0,1,1), Vector3(0,0,-1), Vector3(0,1,0), aspect_ratio, 90, 0.0f, 1.5f);

    std::vector<Sphere> spheres;

    if (!no_output)
        fprintf(stderr, "Output file: %s\n", file_name);
    else {
        fprintf(stderr, "No output will be written\n");
    }

    readInput();
    create_random_scene(spheres);

    auto image_data = static_cast<int*>(malloc(width * height * sizeof(int) * 3));

    // checksums for each color individually
    Checksum checksum(0, 0, 0);
    //pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    // Iterate over each pixel and trace a ray to calculate the color.
    // This is done for samples amount of time for each pixel.
    // TODO: Try to parallelize this.

    //Fai un ciclo che inizializza i thread in base a NUM_THREADS
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    for(int i = 0; i < NUM_THREADS; i++) {
        ThreadData single_data(i, height, width, samples, depth, &camera, &spheres, image_data, &checksum);
        thread_data[i] =single_data;
        pthread_create(&threads[i], NULL, perform_work, &thread_data[i] );
    }
    

    

    // Fai il join prima di salvare

    for (int i = 0; i < NUM_THREADS; i++){
        pthread_join(threads[i], NULL);
    }

    checksum = Checksum(0,0,0); 
    for(int y = 0; y < height; y++) {
            for(int x = 0; x < width; x++) {
                int pos = (y * width + x) * 3;
                checksum += Checksum(image_data[pos],image_data[pos+1],image_data[pos+2]);
            }
        }
    writeOutput(checksum);
    //Saving the render with PPM format
    if(!no_output) {
        FILE* file;
        if ((file = fopen(file_name, "w")) == NULL )
        {
            perror("fopen");
            exit(EXIT_FAILURE);
        }
        if (fprintf(file, "P3\n%d %d %d\n", width, height, 255) < 0)
        {
            perror("fprintf");
            exit(EXIT_FAILURE);
        }
        for(int y = 0; y < height; y++) {
            for(int x = 0; x < width; x++) {
                int pos = (y * width + x) * 3;
                if (fprintf(file, "%d %d %d\n", image_data[pos] , image_data[pos + 1] , image_data[pos + 2] ) < 0)
                {
                    perror("fprintf");
                    exit(EXIT_FAILURE);
                }
            }
        }
        fclose(file);
    }
    
    free(image_data);
    return 0;
}
