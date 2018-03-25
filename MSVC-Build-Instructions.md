## Prepare your Python environment
Run the following commands as administrator.
If your python install is not in your PATH, you need to provide an absolute path to python.exe.
Alternatively, if you are using the Python installed with Visual Studio, you can run dev_env.bat from
the respository, and it will populate PATH with the default path for Python 3.6.
 - python -mpip install -U pip
 - python -mpip install -U matplotlib

## Checking out the source
In a Git Bash window:
 - git clone https://github.com/mclaughc/ppproject.git
 - cd ppproject
 - git submodule init
 - git submodule update

## Compiling the source
 - Open Visual Studio, File->Open->CMake...
 - Locate ppproject\ppproject\native\CMakeLists.txt.
 - Select x64-Release as the configuration, and click Generate.
 - Click CMake->Build All.
 - Click CMake->Install->ppproject_native.

## Running a test program
dev_env.bat assumes 
 - Open a cmd window
 - cd path\to\ppproject
 - dev_env.bat
 - python tests\spectrogrampipeline.py
