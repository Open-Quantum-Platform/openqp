# NMR CGO/GIAO benchmark matrix

Status: reference_ready_openqp_giao_gated

|System|Basis|Method|Backend|Gauge|Status|Max origin delta/ppm|Reference delta/ppm|Wall/s|Peak memory/MiB|Notes|
|---|---:|---|---|---|---|---:|---:|---:|---:|---|
|h2o|sto-3g|hf|pyscf|cgo|ok|343.193360|0.000000|0.140|115.2|PySCF trusted reference|
|h2o|sto-3g|hf|pyscf|giao|ok|0.000253|0.000000|0.115|115.5|PySCF trusted reference|
|h2o|sto-3g|hf|openqp|cgo|ok|not_sampled|pending|0.880|147.4||
|h2o|sto-3g|hf|openqp|giao|gated|pending|pending|0.788|136.3|OpenQP native GIAO path is still NotImplemented|
|h2o|sto-3g|pbe|pyscf|cgo|ok|338.191722|0.000000|0.524|185.4|PySCF trusted reference|
|h2o|sto-3g|pbe|pyscf|giao|ok|9.756080|0.000000|0.516|192.6|PySCF trusted reference|
|h2o|sto-3g|pbe|openqp|cgo|ok|not_sampled|pending|0.811|152.7||
|h2o|sto-3g|pbe|openqp|giao|gated|pending|pending|0.802|152.1|OpenQP native GIAO path is still NotImplemented|
|formaldehyde|6-31g*|hf|pyscf|cgo|ok|108.921635|0.000000|0.818|193.2|PySCF trusted reference|
|formaldehyde|6-31g*|hf|pyscf|giao|ok|0.000004|0.000000|1.054|193.5|PySCF trusted reference|
|formaldehyde|6-31g*|hf|openqp|cgo|ok|not_sampled|pending|0.915|194.5||
|formaldehyde|6-31g*|hf|openqp|giao|gated|pending|pending|0.773|149.4|OpenQP native GIAO path is still NotImplemented|
|formaldehyde|6-31g*|pbe0|pyscf|cgo|ok|104.750288|0.000000|1.781|305.6|PySCF trusted reference|
|formaldehyde|6-31g*|pbe0|pyscf|giao|ok|8.497257|0.000000|2.008|307.4|PySCF trusted reference|
|formaldehyde|6-31g*|pbe0|openqp|cgo|ok|not_sampled|pending|0.822|204.9||
|formaldehyde|6-31g*|pbe0|openqp|giao|gated|pending|pending|0.770|189.8|OpenQP native GIAO path is still NotImplemented|
|formaldehyde|6-31g*|pbe|pyscf|cgo|ok|102.945576|0.000000|1.488|307.5|PySCF trusted reference|
|formaldehyde|6-31g*|pbe|pyscf|giao|ok|10.296883|0.000000|1.728|307.7|PySCF trusted reference|
|formaldehyde|6-31g*|pbe|openqp|cgo|ok|not_sampled|pending|0.786|195.4||
|formaldehyde|6-31g*|pbe|openqp|giao|gated|pending|pending|0.766|193.8|OpenQP native GIAO path is still NotImplemented|
