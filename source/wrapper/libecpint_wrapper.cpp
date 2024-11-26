// libecpint_wrapper.cpp
#include "api.hpp"
#include <vector>
#include <iostream>
#include <memory>
#include <iomanip>

extern "C" {

    typedef struct {
        double *data;
        int size;
    } ResultArray;

    // Initialize integrator and set Gaussian basis
    void* init_integrator(int num_gaussians, double* g_coords, double* g_exps, double* g_coefs, int* g_ams, int* g_lengths) {
        int total_length = 0;
        for (int i = 0; i < num_gaussians; i++) {
            total_length += g_lengths[i];
        }
        printf("Number of Gaussians: %d\n", num_gaussians);


        // Print the values in g_coords
        printf("Gaussian coordinates:\n");
        for (int i = 0; i < num_gaussians * 3; i++) {  // Assuming g_coords has 3 values per Gaussian (x, y, z)
            printf("%f ", g_coords[i]);
            if ((i + 1) % 3 == 0) printf("\n");  // New line after each coordinate triplet
        }

        // Print the values in g_exps
        printf("Gaussian exponents:\n");
        for (int i = 0; i < total_length; i++) {
            printf("%f ", g_exps[i]);
        }
        printf("\n");

        // Print the values in g_coefs
        printf("Gaussian coefficients:\n");
        for (int i = 0; i < total_length; i++) {
            printf("%f ", g_coefs[i]);
        }
        printf("\n");

        // Print the values in g_ams
        printf("Angular momentum numbers:\n");
        for (int i = 0; i < num_gaussians; i++) {
            printf("%d ", g_ams[i]);
        }
        printf("\n");

        // Print the values in g_lengths
        printf("Lengths of Gaussian functions:\n");
        for (int i = 0; i < num_gaussians; i++) {
            printf("%d ", g_lengths[i]);
        }
        printf("\n");
        libecpint::ECPIntegrator* factory = new libecpint::ECPIntegrator();
        factory->set_gaussian_basis(num_gaussians, g_coords, g_exps, g_coefs, g_ams, g_lengths);
        return factory;
    }

    // Set ECP basis
    void set_ecp_basis(void* integrator, int num_ecps, double* u_coords, double* u_exps, double* u_coefs, int* u_ams, int* u_ns, int* u_lengths) {
        libecpint::ECPIntegrator* factory = static_cast<libecpint::ECPIntegrator*>(integrator);
        factory->set_ecp_basis(num_ecps, u_coords, u_exps, u_coefs, u_ams, u_ns, u_lengths);
    }

    // Initialize integrator
    void init_integrator_instance(void* integrator) {
        libecpint::ECPIntegrator* factory = static_cast<libecpint::ECPIntegrator*>(integrator);
        factory->init();
    }

    // Compute integrals and return result in a structure
    ResultArray compute_integrals(void* integrator) {
        libecpint::ECPIntegrator* factory = static_cast<libecpint::ECPIntegrator*>(integrator);
        factory->compute_integrals();
        std::shared_ptr<std::vector<double>> ints = factory->get_integrals();
        ResultArray result;
        result.size = ints->size();
//	printf("%d ", result.size);

        result.data = new double[result.size];
        std::copy(ints->begin(), ints->end(), result.data);
//        std::cout << "result.data contains:\n";
//        for (size_t i = 0; i < 35; ++i) {
//            for (size_t j = 0; j <= i; ++j){
//                std::cout << result.data[i*35+j] << "  \n";
//            }
//        }
//        std::cout << std::endl;
        return result;
    }

    // Clean up
    void free_integrator(void* integrator) {
        delete static_cast<libecpint::ECPIntegrator*>(integrator);
    }

    void free_result(ResultArray result) {
        delete[] result.data;
    }
}

