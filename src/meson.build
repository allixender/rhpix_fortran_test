project('rhpix', 'fortran', meson_version: '>=0.49', version: '1.1.1')

if meson.get_compiler('fortran').get_id() == 'gcc'
  f_args = ['-Wno-maybe-uninitialized']
else
  f_args = []
endif

prefix = get_option('prefix')
libdir = get_option('libdir')

deps = []

rhpix_modules = ['mod_rhpix.f90']

# rhpix_sources = []
# foreach m: rhpix_modules
#   rhpix_sources += [join_paths(src_dir, m)]
# endforeach

rhpix = shared_library('rhpix',
                      rhpix_modules,
                      fortran_args: f_args,
                      dependencies: deps)

rhpix_dep = declare_dependency(link_with : rhpix)

executable(
  'mainApp',
  ['main.f90'],
  fortran_args: f_args,
  link_with: rhpix)

# executable('mainApp', ['main.f90','mod_rhpix.f90'])
