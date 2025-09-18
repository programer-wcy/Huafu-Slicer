A slicer for 3D printing with our own devices. This version of our projects is based on BambuStudio-02.00.03.54

1. unzip BambuStudio_dep_external_src/deps_src.zip to some specific directory (d:\deps_src, for example)

2. Install visual studio 2022 with MSVC of v143
	visualstudiosetup.exe --channelUri https://aka.ms/vs/17/release.LTSC.17.8/channel

3. Change work derectory to ./BambuStudio/deps and excute the following commands
	mkdir build
	mkdir BambuStudio_dep
	cd build
	cmake ../ -G "Visual Studio 17 2022" -A x64 -DDESTDIR="../../../BambuStudio_dep" -DCMAKE_BUILD_TYPE=Release

4. Copy all the unzipped contents described in step 1 to directory ./BambuStudio/deps/build/

5. Open ./BambuStudiodeps/build/BambuStudio-deps.sln with Visual Studio 2022 and build project ALL_BUILD to generate dep libs and headers for the use of later building of BambuStudio
