/* Test-only OpenBLAS control API. Numerical BLAS stays on the native backend.
 * The callback reproduces a setter that also changes the OpenMP width. */
static int width = 4, caps = 0, calls = 0;
static void (*set_omp)(int) = 0;

int openblas_get_num_threads(void) { return width; }
void openblas_set_num_threads(int n)
{
    width = n;
    calls++;
    if (n == 1) caps++;
    if (set_omp) set_omp(n);
}
void ecp_fixture_callback(void (*fn)(int)) { set_omp = fn; }
int ecp_fixture_caps(void) { return caps; }
int ecp_fixture_calls(void) { return calls; }
