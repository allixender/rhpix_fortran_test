project('rhpix',  ['fortran','c'], meson_version: '>=0.49', version: '1.1.1')

if meson.get_compiler('fortran').get_id() == 'gcc'
  f_args = ['-Wno-maybe-uninitialized']
else
  f_args = []
endif

deps = []

rhpix_modules = ['mod_rhpix_utils.f90']

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

executable('test1c',
    ['test1.c'],
    c_args: ['-fPIC', '-lrhpix', '-L.'],
    link_args: ['-lrhpix', '-L.'],
    dependencies: rhpix_dep)
