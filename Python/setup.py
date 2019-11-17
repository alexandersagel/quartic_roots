from distutils.core import setup, Extension

module1 = Extension('quartic_c',
                    sources = ['quarticRoots.cpp'])

setup (name = 'Quartic',
       version = '1.0',
       description = 'Computing zeros of quartic polynomials',
       ext_modules = [module1])