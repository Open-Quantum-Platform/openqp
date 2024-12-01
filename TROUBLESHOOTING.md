# Compilation

## MacOS compilation

### ld: Assertion failed: (resultIndex < sectData.atoms.size()), function findAtom, file Relocations.cpp, line 1336.

Depending on MacOSX SDK, compilation may stop at the linking stage with the following message:
```txt
0  0x10e045f63  __assert_rtn + 64
1  0x10df47f63  ld::AtomPlacement::findAtom(unsigned char, unsigned long long, ld::AtomPlacement::AtomLoc const*&, long long&) const + 1411
2  0x10df64451  ld::InputFiles::SliceParser::parseObjectFile(mach_o::Header const*) const + 19745
3  0x10df7190a  ld::InputFiles::SliceParser::parse() const + 3242
4  0x10df74b91  ld::InputFiles::parseAllFiles(void (ld::AtomFile const*) block_pointer)::$_7::operator()(unsigned long, ld::FileInfo const&) const + 657
5  0x7ff804cd0066  _dispatch_client_callout2 + 8
6  0x7ff804ce1e09  _dispatch_apply_invoke + 213
7  0x7ff804cd0033  _dispatch_client_callout + 8
8  0x7ff804ce00f6  _dispatch_root_queue_drain + 683
9  0x7ff804ce0768  _dispatch_worker_thread2 + 170
10  0x7ff804e6dc0f  _pthread_wqthread + 257
ld: Assertion failed: (resultIndex < sectData.atoms.size()), function findAtom, file Relocations.cpp, line 1336.
```

To resolve this issue, one needs to pass `-Wl,-ld_classic` flags to linker. With CMake, extra flags to configure command should be passed:
```bash
-DCMAKE_C_FLAGS="-Wl,-ld_classic" -DCMAKE_CXX_FLAGS="-Wl,-ld_classic" -DCMAKE_Fortran_FLAGS="-Wl,-ld_classic"
```

If you use additional flags in your configuration command, adjust it in a proper way.
