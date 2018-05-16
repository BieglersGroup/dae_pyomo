from setuptools import setup

setup(
    name='dae_pyomo',
    version='1.0',
    packages=['capd_short', 'capd_short.gams', 'capd_short.dae_v2', 'capd_short.dae_v2.collocation_funcs',
              'capd_short.pyomodae'],
    url='https://github.com/BieglersGroup',
    license='MIT',
    author='David Thierry',
    author_email='',
    description='Examples of modelling with DAEs'
)
