import os
import sys
import subprocess
from pathlib import Path
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cfg = "DEBUG" if self.debug else "RELEASE"

        build_args = []
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(extdir),
            "-DCMAKE_BUILD_TYPE={}".format(cfg)
        ]
        cmake_args += ["-DPYTHON=ON"]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )
        # Move from build temp to final position
        for ext in self.extensions:
            self.move_output(ext)

    def move_output(self, ext):
        build_temp = Path(self.build_temp).resolve()
        dest_path = Path(self.get_ext_fullpath(ext.name)).resolve()
        source_path = build_temp / "lib" / self.get_ext_filename(ext.name)
        dest_directory = dest_path.parents[0]
        dest_directory.mkdir(parents=True, exist_ok=True)
        self.copy_file(source_path, dest_path)

if not os.path.isfile("src/pyommp_interface.cpp"):
    raise RuntimeError("Running setup.py is only supported "
                       "from top level of repository as './setup.py <command>'")
def ommp_v(version):
    return str(version.tag)

def ommp_local_v(version):
    if version.distance > 0 or version.dirty:
        lv = '+'
    else:
        lv = ''
    if version.distance > 0:
        lv += "r{:d}.{:s}".format(version.distance, version.node)
        if version.dirty:
            lv += '.dirty'
    elif version.dirty:
         lv += 'dirty'
    return lv

setup(
    name="pyopenmmpol",
    ext_modules=[CMakeExtension("pyopenmmpol")],
    cmdclass={"build_ext": CMakeBuild},
    use_scm_version = {
        "root": "./",
        "relative_to": __file__,
        "version_scheme": ommp_v,
        "local_scheme": ommp_local_v,
    },
)
