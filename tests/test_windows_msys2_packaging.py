import unittest
import re
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class WindowsMsys2PackagingTests(unittest.TestCase):
    def test_windows_packaging_docs_drop_gui_and_gpu_scope(self):
        docs = (ROOT / "docs" / "windows-msys2-packaging.md").read_text()

        self.assertIn("without GUI", docs)
        self.assertIn("No GPU or CUDA support", docs)
        self.assertIn("MSYS2 UCRT64", docs)
        self.assertIn("LINALG_LIB=netlib", docs)
        self.assertIn("PKGBUILD", docs)
        self.assertIn("pacman -U", docs)
        self.assertIn("basis_set_exchange", docs)
        self.assertIn("GitHub Release", docs)
        self.assertIn("release: published", docs)
        self.assertIn("double-click", docs)
        self.assertIn("install_openqp_msys2.cmd", docs)
        self.assertIn("OpenQP-MSYS2-UCRT64-<version>.msi", docs)
        self.assertIn("preview artifact", docs)
        self.assertIn("clean-machine", docs)
        self.assertIn("pacman -R mingw-w64-ucrt-x86_64-openqp", docs)
        self.assertIn("WiX", (ROOT / "tools" / "windows_msys2" / "README.md").read_text())

    def test_msys2_source_script_uses_ucrt64_and_cli_validation(self):
        script = (ROOT / "tools" / "windows_msys2" / "build_msys2.sh").read_text()

        self.assertIn('MSYSTEM:-}" != "UCRT64"', script)
        self.assertIn("git \\", script)
        self.assertIn("mingw-w64-ucrt-x86_64-gcc-fortran", script)
        self.assertIn("-DLINALG_LIB=netlib", script)
        self.assertIn('pip install --user --break-system-packages --upgrade', script)
        self.assertIn('"jsonschema<4.18"', script)
        self.assertIn('external_root="$(cygpath -m "$repo_root/.cache/openqp/externals")"', script)
        self.assertIn("python -m build --wheel --no-isolation --skip-dependency-check", script)
        self.assertIn("PYTHONPATH", script)
        self.assertIn('$prefix/bin:$PATH', script)
        self.assertIn("activate-openqp.sh", script)
        self.assertIn("liboqp.dll", script)
        self.assertIn('openqp "$example" --omp 2', script)

    def test_pkgbuild_defines_real_msys2_package(self):
        pkgbuild = (
            ROOT / "tools" / "windows_msys2" / "mingw-w64-openqp" / "PKGBUILD"
        ).read_text()
        install = (
            ROOT / "tools" / "windows_msys2" / "mingw-w64-openqp" / "openqp.install"
        ).read_text()
        builder = (ROOT / "tools" / "windows_msys2" / "build_pkg.sh").read_text()

        self.assertIn('pkgbase=mingw-w64-${_realname}', pkgbuild)
        self.assertIn('pkgname=("${MINGW_PACKAGE_PREFIX}-${_realname}")', pkgbuild)
        self.assertIn("pyproject.toml", pkgbuild)
        self.assertIn("pkgver=$(", pkgbuild)
        self.assertIn('mingw_arch=("ucrt64")', pkgbuild)
        self.assertIn('${MINGW_PACKAGE_PREFIX}-python-installer', pkgbuild)
        self.assertIn('_external_root="$(cygpath -m "${srcdir}/openqp-externals")"', pkgbuild)
        self.assertIn('python -m build --wheel --no-isolation --skip-dependency-check', pkgbuild)
        self.assertIn('python -m installer --destdir="${pkgdir}"', pkgbuild)
        self.assertIn("basis_set_exchange", install)
        self.assertIn("makepkg-mingw", builder)
        self.assertIn("sed -i 's/\\r$//' PKGBUILD ./*.install", builder)
        self.assertIn("pacman -U --noconfirm", builder)
        self.assertIn('--break-system-packages --upgrade "jsonschema<4.18" basis_set_exchange geometric', builder)

    def test_clickable_installer_wraps_real_pacman_package(self):
        cmd = (ROOT / "tools" / "windows_msys2" / "install_openqp_msys2.cmd").read_text()
        ps1 = (ROOT / "tools" / "windows_msys2" / "install_openqp_msys2.ps1").read_text()

        self.assertIn("powershell.exe -NoProfile -ExecutionPolicy Bypass", cmd)
        self.assertIn("install_openqp_msys2.ps1", cmd)
        self.assertIn("pause", cmd.lower())
        self.assertIn("Find-Msys2Root", ps1)
        self.assertIn("C:\\msys64", ps1)
        self.assertIn("usr\\bin\\pacman.exe", ps1)
        self.assertNotIn("ucrt64\\bin\\python.exe", ps1)
        self.assertIn("mingw-w64-ucrt-x86_64-python", ps1)
        self.assertIn("mingw-w64-ucrt-x86_64-python-cffi", ps1)
        self.assertIn("mingw-w64-ucrt-x86_64-python-numpy", ps1)
        self.assertIn("mingw-w64-ucrt-x86_64-python-scipy", ps1)
        self.assertIn("mingw-w64-ucrt-x86_64-openqp-*.pkg.tar.*", ps1)
        self.assertIn("pacman -U --noconfirm", ps1)
        self.assertIn('--break-system-packages --upgrade "jsonschema<4.18" basis_set_exchange geometric', ps1)
        self.assertIn("OpenQP UCRT64.cmd", ps1)
        self.assertIn("openqp-msys2-ucrt64.pkg.tar.*", ps1)
        self.assertIn("command -v openqp", ps1)

    def test_msi_recipe_bundles_click_installer_and_pacman_package(self):
        wxs_path = ROOT / "tools" / "windows_msys2" / "msi" / "Product.wxs"
        build_msi = (ROOT / "tools" / "windows_msys2" / "build_msi.ps1").read_text()
        wxs = wxs_path.read_text()

        ET.parse(wxs_path)
        self.assertIn("OpenQP MSYS2 UCRT64 Installer", wxs)
        self.assertIn("UpgradeCode=", wxs)
        self.assertIn("install_openqp_msys2.cmd", wxs)
        self.assertIn("install_openqp_msys2.ps1", wxs)
        self.assertIn("openqp-msys2-ucrt64.pkg.tar.zst", wxs)
        self.assertIn("Install OpenQP for MSYS2", wxs)
        self.assertIn("dotnet tool install --global wix", (ROOT / ".github" / "workflows" / "windows-msys2.yml").read_text())
        self.assertIn("OpenQP-MSYS2-UCRT64-$ProductVersion.msi", build_msi)
        self.assertIn("pyproject.toml", build_msi)
        self.assertIn("wix", build_msi)

    def test_pkgbuild_version_matches_project_version(self):
        pyproject = (ROOT / "pyproject.toml").read_text()
        pkgbuild = (
            ROOT / "tools" / "windows_msys2" / "mingw-w64-openqp" / "PKGBUILD"
        ).read_text()

        project_version = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, re.M).group(1)
        command = re.search(r"^pkgver=\$\((.+)\)", pkgbuild, re.M).group(1)
        self.assertIn("pyproject.toml", command)
        self.assertIn(project_version, pyproject)

    def test_manual_workflow_is_separate_from_release_wheel_matrix(self):
        workflow = (ROOT / ".github" / "workflows" / "windows-msys2.yml").read_text()
        release_workflow = (ROOT / ".github" / "workflows" / "build_wheels.yml").read_text()

        self.assertIn("workflow_dispatch", workflow)
        self.assertIn("release:", workflow)
        self.assertIn("types: [published]", workflow)
        self.assertIn("contents: write", workflow)
        self.assertIn("msys2/setup-msys2", workflow)
        self.assertIn("msystem: UCRT64", workflow)
        self.assertIn("build_pkg.sh", workflow)
        self.assertIn("install_openqp_msys2.cmd", workflow)
        self.assertIn("install_openqp_msys2.ps1", workflow)
        self.assertIn("build_msi.ps1", workflow)
        self.assertIn("Validate MSI wrapper install", workflow)
        self.assertIn("msiexec.exe", workflow)
        self.assertIn("OpenQP MSYS2 Installer", workflow)
        self.assertIn("*.msi", workflow)
        self.assertIn("*.pkg.tar.*", workflow)
        self.assertIn("SHA256SUMS.txt", workflow)
        self.assertIn("gh release upload", workflow)
        self.assertIn("github.event.release.tag_name", workflow)
        self.assertNotIn("windows-msys2", release_workflow)
        self.assertNotIn("msys2/setup-msys2", release_workflow)


if __name__ == "__main__":
    unittest.main()
