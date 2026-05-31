# NMR CGO/GIAO benchmark matrix

Status: reference_ready_openqp_giao_gated

|System|Basis|Method|Backend|Gauge|Status|Max origin delta/ppm|Reference delta/ppm|Wall/s|Peak memory/MiB|Notes|
|---|---:|---|---|---|---|---:|---:|---:|---:|---|
|h2o|sto-3g|hf|pyscf|cgo|ok|343.193360|0.000000|0.171|141.9|PySCF trusted reference|
|h2o|sto-3g|hf|pyscf|giao|ok|0.000253|0.000000|0.115|142.1|PySCF trusted reference|
|h2o|sto-3g|hf|openqp|cgo|ok|not_sampled|pending|0.837|143.1||
|h2o|sto-3g|hf|openqp|giao|gated|pending|pending|0.764|136.3|OpenQP native GIAO path is still NotImplemented|
|h2o|sto-3g|pbe|pyscf|cgo|ok|338.191722|0.000000|0.526|213.4|PySCF trusted reference|
|h2o|sto-3g|pbe|pyscf|giao|ok|9.756080|0.000000|0.514|215.7|PySCF trusted reference|
|h2o|sto-3g|pbe|openqp|cgo|ok|not_sampled|pending|0.766|155.0||
|h2o|sto-3g|pbe|openqp|giao|gated|pending|pending|0.762|147.2|OpenQP native GIAO path is still NotImplemented|
|formaldehyde|6-31g*|hf|pyscf|cgo|ok|108.921635|0.000000|0.838|216.3|PySCF trusted reference|
|formaldehyde|6-31g*|hf|pyscf|giao|ok|0.000004|0.000000|1.068|216.4|PySCF trusted reference|
|formaldehyde|6-31g*|hf|openqp|cgo|ok|not_sampled|pending|0.915|195.5||
|formaldehyde|6-31g*|hf|openqp|giao|gated|pending|pending|0.763|148.3|OpenQP native GIAO path is still NotImplemented|
|formaldehyde|6-31g*|pbe0|pyscf|cgo|ok|104.750288|0.000000|1.772|332.5|PySCF trusted reference|
|formaldehyde|6-31g*|pbe0|pyscf|giao|ok|8.497257|0.000000|2.017|332.6|PySCF trusted reference|
|formaldehyde|6-31g*|pbe0|openqp|cgo|ok|not_sampled|pending|0.811|210.7||
|formaldehyde|6-31g*|pbe0|openqp|giao|gated|pending|pending|0.776|189.3|OpenQP native GIAO path is still NotImplemented|
|formaldehyde|6-31g*|pbe|pyscf|cgo|ok|102.945576|0.000000|1.529|332.8|PySCF trusted reference|
|formaldehyde|6-31g*|pbe|pyscf|giao|ok|10.296883|0.000000|1.779|336.0|PySCF trusted reference|
|formaldehyde|6-31g*|pbe|openqp|cgo|ok|not_sampled|pending|0.808|201.5||
|formaldehyde|6-31g*|pbe|openqp|giao|gated|pending|pending|0.805|192.4|OpenQP native GIAO path is still NotImplemented|
