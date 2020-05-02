from setuptools import setup

setup(
   name='HFSSdrawpy',
   version='0.1',
   description='A module for drawing circuits in HFSS with python',
   author='HQC',
   # author_email='foomail@foo.com',
   packages=['HFSSdrawpy'],  #same as name
   install_requires=['bar', 'greek'], #external packages as dependencies
   # scripts=[
   #          'scripts/cool',
   #          'scripts/skype',
   #         ]
)