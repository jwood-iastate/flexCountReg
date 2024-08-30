Note:
To compile the c++ code, the GNU GSL library is required. To install the GSL library:

On Windows, you can install it from (assuming yuo have Visual Studio with c++ development installed):
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
.\vcpkg integrate install
.\vcpkg install gsl


Intall on Linux: 
sudo apt-get install libgsl-dev


Install on Mac: 
brew install gsl