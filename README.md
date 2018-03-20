# Modern C++ Continuous Integration

## Pre-installation

```shell
  sudo apt install catch
```


## Generate project

```shell
  cmake -H. -BBuild
```

Auto detect everything.

If you like to set a implicit compiler set the variable CXX=${COMPILER}, for example COMPILER could be gcc, clang and so on.

Auto detect in Windows usually generate a Visual Studio project since msbuild require it, but in OSX does not generate and XCode project, since is not required for compiling using XCode clang.

Specify build type debug/release

```shell
  # generate a debug project
  cmake -H. -BBuild -DCMAKE_BUILD_TYPE=Debug
  # generate a release project
  cmake -H. -BBuild -DCMAKE_BUILD_TYPE=Release
```

Specify architecture

```shell
  # 64 bits architecture
  cmake -H. -BBuild -Ax64
  # ARM architecture
  cmake -H. -BBuild -AARM
  # Windows 32 bits architecture
  cmake -H. -BBuild -AxWin32
```

Generate different project types

```shell
  # MinGW makefiles
  cmake -H. -BBuild -G "MinGW Makefiles"
  # XCode project
  cmake -H. -BBuild -G "XCode"
  # Visual Studio 15 2017 solution
  cmake -H. -BBuild -G "Visual Studio 15 2017"
```

## Build

From the Build folder

```shell
  # build the default build type (in multi build types usually debug)
  cmake --build .
  # build a specific build type
  cmake --build . --config Release
```
## Run tests

From the Build folder

```shell
  # run all test using the default build type
  ctest -V
  # run all test in Release build type
  ctest -V -C Release
```

## references

- https://cmake.org/
- https://docs.travis-ci.com/user/languages/cpp/
- https://www.appveyor.com/docs/lang/cpp/
- https://github.com/isocpp/CppCoreGuidelines
- https://github.com/philsquared/Catch
- https://github.com/gabime/spdlog
- https://github.com/cognitivewaves/CMake-VisualStudio-Example
- http://derekmolloy.ie/hello-world-introductions-to-cmake/
- https://cmake.org/Wiki/CMake/Testing_With_CTest
- https://www.appveyor.com/docs/lang/cpp/
- https://docs.travis-ci.com/user/languages/cpp/
- https://github.com/philsquared/Catch/blob/master/docs/build-systems.md
- https://stackoverflow.com/questions/14446495/cmake-project-structure-with-unit-tests
