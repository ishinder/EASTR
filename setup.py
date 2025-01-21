"""EASTR build."""

import os
import platform
import subprocess

from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
  def __init__(self, name, sourcedir=''):
    Extension.__init__(self, name, sources=[])
    self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
  """CMake build extension."""
  def run(self):
    extensions = ', '.join(e.name for e in self.extensions)
    deps = ['cmake', 'make']
    for dep in deps:
      try:
        subprocess.check_output([dep, '--version'])
      except OSError as e:
        raise RuntimeError(
          f'{dep} must be installed to build the following extensions: {extensions}') from e

    for ext in self.extensions:
      self.build_extension(ext)

  def build_extension(self, ext):
    extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
    extdir = os.path.join(extdir, ext.name)
    cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir]
    cmake_args += ['-DEXECUTABLE_OUTPUT_PATH=' + extdir]

    cfg = 'Debug' if self.debug else 'Release'
    build_args = ['--config', cfg]

    if platform.system() == 'Windows':
      cmake_args += ['-DCMAKE_GENERATOR_PLATFORM=x64']
      cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
      build_args += ['--', '/m']
    else:
      cmake_args += ['-DCMAKE_BUILD_TYPE=Release .' + cfg]

    env = os.environ.copy()
    env['CXXFLAGS'] = '{} -DVERSION_INFO=\"{}\"'.format(
      env.get('CXXFLAGS', ''),
      self.distribution.get_version())
    if not os.path.exists(self.build_temp):
      os.makedirs(self.build_temp)
    subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                          cwd=self.build_temp, env=env)
    subprocess.check_call(['cmake', '--build', '.', *build_args],
                          cwd=self.build_temp)

desc = 'Tool for emending alignments of spuriously spliced transcript reads'

with open('./README.md', 'r', encoding='utf-8') as fh:
  long_description = fh.read()

with open('./LICENSE', 'r', encoding='utf-8') as fh:
  license_str = fh.read()

setup(
    name='eastr',
    version='1.1.1',  # https://semver.org/
    author='Ida Shinder',
    author_email='ishinde1@jhmi.edu',
    description=desc,
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/ishinder/EASTR',
    cmdclass=dict(build_ext=CMakeBuild),
    ext_modules=[CMakeExtension('eastr', 'utils')],
    install_requires=[
      'biopython>=1.81,<2.0',
      'mappy>=2.26,<3.0',
      'numpy>=1.26.1',
      'pandas>=2.1.2,<2.3',
      'pysam>=0.22.0,<0.23',
    ],
    packages=find_packages(
      where='./',
      include=['EASTR'],
    ),
    entry_points={
      'console_scripts': [
          'eastr = EASTR.run_eastr:main',
      ]
    },
    python_requires='>=3.10',
    license='MIT',
    license_files=('LICENSE',),
)
