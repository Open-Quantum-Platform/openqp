# NMR CGO/GIAO benchmark matrix

Status: reference_ready_openqp_giao_gated

|System|Basis|Method|Backend|Gauge|Status|Max origin delta/ppm|Reference delta/ppm|Wall/s|Notes|
|---|---:|---|---|---|---|---:|---:|---:|---|
|h2o|sto-3g|hf|pyscf|cgo|ok|343.193360|0.000000|0.140|PySCF trusted reference|
|h2o|sto-3g|hf|pyscf|giao|ok|0.000253|0.000000|0.114|PySCF trusted reference|
|h2o|sto-3g|hf|openqp|cgo|skipped|not_sampled|pending|0.000|OPENQP_ROOT does not point to a built OpenQP package|
|h2o|sto-3g|hf|openqp|giao|skipped|pending|pending|0.000|OPENQP_ROOT does not point to a built OpenQP package|
|h2o|sto-3g|pbe|pyscf|cgo|ok|338.191722|0.000000|0.514|PySCF trusted reference|
|h2o|sto-3g|pbe|pyscf|giao|ok|9.756080|0.000000|0.519|PySCF trusted reference|
|h2o|sto-3g|pbe|openqp|cgo|skipped|not_sampled|pending|0.000|OPENQP_ROOT does not point to a built OpenQP package|
|h2o|sto-3g|pbe|openqp|giao|skipped|pending|pending|0.000|OPENQP_ROOT does not point to a built OpenQP package|
|formaldehyde|6-31g*|hf|pyscf|cgo|ok|108.921635|0.000000|0.805|PySCF trusted reference|
|formaldehyde|6-31g*|hf|pyscf|giao|ok|0.000004|0.000000|1.080|PySCF trusted reference|
|formaldehyde|6-31g*|hf|openqp|cgo|skipped|not_sampled|pending|0.000|OPENQP_ROOT does not point to a built OpenQP package|
|formaldehyde|6-31g*|hf|openqp|giao|skipped|pending|pending|0.000|OPENQP_ROOT does not point to a built OpenQP package|
|formaldehyde|6-31g*|pbe0|pyscf|cgo|ok|104.750288|0.000000|1.801|PySCF trusted reference|
|formaldehyde|6-31g*|pbe0|pyscf|giao|ok|8.497257|0.000000|2.045|PySCF trusted reference|
|formaldehyde|6-31g*|pbe0|openqp|cgo|skipped|not_sampled|pending|0.000|OPENQP_ROOT does not point to a built OpenQP package|
|formaldehyde|6-31g*|pbe0|openqp|giao|skipped|pending|pending|0.000|OPENQP_ROOT does not point to a built OpenQP package|
|formaldehyde|6-31g*|pbe|pyscf|cgo|ok|102.945576|0.000000|1.513|PySCF trusted reference|
|formaldehyde|6-31g*|pbe|pyscf|giao|ok|10.296883|0.000000|1.798|PySCF trusted reference|
|formaldehyde|6-31g*|pbe|openqp|cgo|skipped|not_sampled|pending|0.000|OPENQP_ROOT does not point to a built OpenQP package|
|formaldehyde|6-31g*|pbe|openqp|giao|skipped|pending|pending|0.000|OPENQP_ROOT does not point to a built OpenQP package|
