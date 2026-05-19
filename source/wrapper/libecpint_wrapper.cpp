/**
 * @file libecpint_c_api.cpp
 * @brief C ABI wrapper around libecpint::ECPIntegrator for use from Fortran/C.
 * @author Mohsen Mazaherifar
 * @date January 2025
 */
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
    void init_integrator_instance(void* integrator, int deriv_order) {
        libecpint::ECPIntegrator* factory = static_cast<libecpint::ECPIntegrator*>(integrator);
        factory->init(deriv_order);
    }

    // Compute integrals and return result in a structure
    ResultArray compute_integrals(void* integrator) {
        libecpint::ECPIntegrator* factory = static_cast<libecpint::ECPIntegrator*>(integrator);
        factory->compute_integrals();
        std::shared_ptr<std::vector<double>> ints = factory->get_integrals();
        ResultArray result;
        result.size = ints->size();

        result.data = new double[result.size];
        std::copy(ints->begin(), ints->end(), result.data);
        return result;
    }

    //Compute first derivs
    ResultArray compute_first_derivs(void* integrator) {
        libecpint::ECPIntegrator* factory = static_cast<libecpint::ECPIntegrator*>(integrator);
        factory->compute_first_derivs();
        std::vector<std::shared_ptr<std::vector<double>>> first_derivs = factory->get_first_derivs();

        int total_size = 0;
        for (const auto& vec : first_derivs) {
            total_size += vec->size();
        }

        ResultArray result;
        result.size = total_size;
        result.data = new double[total_size];
        int index = 0;
        for (const auto& vec : first_derivs) {
            std::copy(vec->begin(), vec->end(), result.data + index);
            index += vec->size();
        }

        return result;
    }

    // Clean up
    void free_integrator(void* integrator) {
        delete static_cast<libecpint::ECPIntegrator*>(integrator);
	integrator = nullptr;
    }

    void free_result(ResultArray result) {
        delete[] result.data;
	result.data = nullptr;
    }
}

