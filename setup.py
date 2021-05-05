from setuptools import setup, find_packages

def readme():
  with open('README.rst', 'r') as fh:
    return fh.read()

setup(
  entry_points={
      'console_scripts': [
          'bernoulli = polylib.bernoulli:cli',
      ],
      'gui_scripts': [],
  },
  name='polylib',
  url='https://github.com/sj-simmons/polylib/archive/v0.1.tar.gz',
  author='Scott Simmons',
  author_email='ssimmons@drury.edu',
  packages=find_packages(),
  python_requires='>=3.6',
  install_requires=[],
  version='0.1',
  license='Apache 2.0',
  description='polynomials over a rings or a field',
  long_description=readme(),
  include_package_data=True,
  zip_safe=False,
  project_urls={'Upstream Repository': 'https://gihub.com/sj-simmons/polylib'},
  classifiers=[
      'Development Status :: 4 - Beta',
      'Intended Audience :: Science/Research',
      'Intended Audience :: Education',
      'Programming Language :: Python :: 3.6',
      'Programming Language :: Python :: 3.7',
      'Programming Language :: Python :: 3.8',
      'Programming Language :: Python :: 3.9',
      'Topic :: Scientific/Engineering',
      'Topic :: Scientific/Engineering :: Mathematics',
      'License :: OSI Approved :: Apache Software License'
  ]
)
