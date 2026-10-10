# Data Tags and Tagarray

OpenQP uses a tagarray container to move arrays and scalar records across
Fortran, C, and Python. It is the shared storage for matrices, orbital data,
response vectors, SOC arrays, Hessians, and other intermediate or result data
that must survive across language boundaries.

The container is provided by
[Open-Quantum-Platform/tagarray](https://github.com/Open-Quantum-Platform/tagarray)
(v1.0.0), pinned as an external dependency in `external/CMakeLists.txt`.

## Naming

Tag names use the `OQP::` prefix. Fortran code usually imports named constants
from `source/tagarray_driver.F90`, such as:

| Tag | Typical content |
| --- | --- |
| `OQP::DM_A`, `OQP::DM_B` | Alpha and beta density matrices. |
| `OQP::FOCK_A`, `OQP::FOCK_B` | Alpha and beta Fock matrices. |
| `OQP::E_MO_A`, `OQP::E_MO_B` | Molecular-orbital energies. |
| `OQP::VEC_MO_A`, `OQP::VEC_MO_B` | Molecular-orbital coefficient matrices. |
| `OQP::Hcore`, `OQP::SM`, `OQP::TM` | Core Hamiltonian, overlap, and kinetic-energy matrices. |
| `OQP::td_energies`, `OQP::td_bvec_mo` | TDHF/TDDFT/MRSF response data. |
| `OQP::soc_eval`, `OQP::soc_hsoc_re`, `OQP::soc_hsoc_im` | SOC eigenvalues and Hamiltonian components. |
| `OQP::hf_hessian` | Native analytical HF/DFT Hessian matrix. |

## Fortran Pattern

### Creating an output record

To create a container-owned array and obtain a writable, correctly typed Fortran
pointer to it, use the one-call `alloc_or_die`:

```fortran
use oqp_tagarray_driver, only: OQP_hf_hessian, OQP_hf_hessian_comment

real(kind=dp), contiguous, pointer :: hess_store(:,:)

! Allocate the record AND bind the pointer in a single call. The TA_TYPE_* id is
! inferred from the pointer kind, the shape is given as plain integers, an
! existing record under the tag is replaced, and a failure aborts with a message.
call infos%dat%alloc_or_die(OQP_hf_hessian, (/ ncart, ncart /), hess_store, &
                            description=OQP_hf_hessian_comment)
```

`alloc_or_die` replaces the older three-step ritual (`reserve_data`, then
`data_has_tags`, then `tagarray_get_data`). A status-returning `alloc` variant is
available when a module wants to handle the error itself, and a free-function
spelling `ta_allocate(infos%dat, ...)` exists for those who prefer it:

```fortran
integer(c_int32_t) :: status
status = infos%dat%alloc(OQP_hf_hessian, (/ ncart, ncart /), hess_store, &
                         description=OQP_hf_hessian_comment)
if (status /= TA_OK) call show_message(get_status_message(status), WITH_ABORT)
```

Supported element types are real64/real32/int32/int64/complex(real64), pointer
ranks 1–3.

### Reading an existing record

To consume a record produced by an earlier stage, confirm it exists and bind a
typed pointer with `tagarray_get_data`:

```fortran
use oqp_tagarray_driver, only: OQP_SM, data_has_tags, tagarray_get_data

real(kind=dp), contiguous, pointer :: smat(:)

call data_has_tags(infos%dat, [character(len=80) :: OQP_SM], &
                   module_name, subroutine_name, abort=.true.)
call tagarray_get_data(infos%dat, OQP_SM, smat)
```

### Low-level operations

For records that `alloc_or_die` does not cover — an `integer(8)`-sized shape
(for example the `nbf**4` in-core ERI tensor) or a tag that is filled on the C
side with no Fortran pointer — use the low-level container methods directly. The
`create` shape is `integer(c_int64_t)`:

```fortran
status = infos%dat%create(OQP_ERI_AO, TA_TYPE_REAL64, [nbf4], &
                          description=OQP_ERI_AO_comment, override=.true.)
call infos%dat%erase([character(len=80) :: OQP_QMAT])   ! invalidate a stale tag
```

Use the typed constants and comments in `tagarray_driver.F90` instead of
hard-coded strings when adding stable records.

## C Bridge

The native C interface exposes generic tag access through:

| Function | Purpose |
| --- | --- |
| `oqp_get(handle, code, type_id, ndims, dims, value)` | Return a pointer, type, dimensions, and data length for an existing tag. |
| `oqp_alloc(handle, code, type_id, ndims, dims, value)` | Allocate storage for a tag and return a writable pointer. |
| `oqp_del(handle, code)` | Delete a tag record. |

These are low-level functions. Callers must respect the returned type and
dimensions.

## Python Access

Python accesses tagarray records through `Molecule.data`:

```python
td_energies = mol.data["OQP::td_energies"]
soc_eval = mol.data["OQP::soc_eval"]
```

`Molecule` also creates convenience getters for common tags. For example,
`OQP::DM_A` becomes `mol.get_dm_a()`, and `OQP::td_energies` becomes
`mol.get_td_energies()`.

## Contributor Checklist

When adding a new shared data record:

1. Add a stable `OQP::...` constant and comment in `source/tagarray_driver.F90`.
2. Create the record and bind its pointer in one call with
   `infos%dat%alloc_or_die(tag, shape, ptr, description=...)` (or the low-level
   `create` for `integer(8)` shapes or C-filled tags), passing the shape that
   matches the pointer rank.
3. Read upstream records with `data_has_tags` + `tagarray_get_data`.
4. Add it to Python result handling only if users or workflows need it.
5. Add tests that check the tag is populated with the expected shape and data
   type.
6. Document user-facing result tags in [Results and Molecule Data](../api/results.md).
