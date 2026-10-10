# XC response scaffold

Preserved from OpenQP `beb51a2617ad`; not included in CMake production targets.
The CUDA operation is a packed elementwise density/kernel product. Full XC
quadrature, derivative channels, and cache ownership are not implemented here.
It must not replace a validated CPU XC response. Planning contracts are in
`python/openqp_gpu/gpu.py` and `tdhf_xc_response_cache.py`.
