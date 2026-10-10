// Declarations for the embedded RTC2 rotation tables (see
// routec_tables_embedded.cpp). Lets grad.cu load the tables from memory via
// fmemopen, so the shipped library needs no external routec_tables.bin.
#pragma once
#include <cstddef>
extern const unsigned char routec_tables_embedded[];
extern const size_t routec_tables_embedded_len;
