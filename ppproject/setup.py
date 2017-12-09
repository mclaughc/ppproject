from distutils.core import setup, Extension
from glob import glob

shared_sources = glob('native/shared/*.cpp')

_splitter = Extension('_splitter',
                      define_macros = [('MACRO', 'VALUE')],
                      include_dirs = ['native'],
                      libraries = ['sndfile'],
                      library_dirs = ['/usr/local/lib'],
                      sources = shared_sources + ['native/splitter/splitter.cpp'])

setup (name = 'PackageName',
       version = '1.0',
       description = 'This is a demo package',
       author = 'Author',
       author_email = 'email@example.com',
       url = 'https://google.com/',
       long_description = '''
Demo package
''',
       ext_modules = [_splitter])


