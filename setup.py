from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='mtunneling',
      version='0.1',
      description='macroscopic tunneling with coherent states',
      long_description=readme(),
      author='Nicolas Gheeraert',
      author_email='n.gheeraert@physics.iitm.ac.in',
      license='',
      packages=['mtunneling'],
      install_requires=[],
      include_package_data=True,
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose'],)
