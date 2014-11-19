from distutils.core import setup, Extension
setup(name="PositionWeightMatrix", version="1.0",
	ext_modules = [
		Extension("PositionWeightMatrix", ["PositionWeightMatrix.c"])
		])
