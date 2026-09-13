# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

from distutils.core import setup, Extension

include_dirs = []
include_dirs.append("/dasch/Pipeline")
include_dirs.append("/dasch/install/include")
include_dirs.append("/usr/include/mysql")
define_macros = []
link_args = []
link_args.append("/dasch/Pipeline/pipelineutils.a")
link_args.append("/dasch/install/lib/libwcs.a")
link_args.append("/dasch/install/lib/libtable.a")
link_args.append("/dasch/install/lib/libutil.a")
lib_dirs = []
lib_dirs.append("/dasch/install/lib")
lib_dirs.append("/usr/lib64/mysql")
extra_link_libs = []
extra_link_libs.append("mysqlclient")

setup(
    name="dasch",
    ext_modules=[
        Extension(
            "dasch",
            include_dirs=include_dirs,
            define_macros=define_macros,
            library_dirs=lib_dirs,
            libraries=extra_link_libs,
            extra_link_args=link_args,
            sources=["pythondasch.c"],
        )
    ],
)
