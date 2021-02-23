
currently (Feb 2021) use one before latest meson to avoid "not mentioned in its dyndep file" error:

pip install ninja meson==0.56

set PATH=C:\dev\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\bin;%PATH%

mkdir build2
cd build2
meson ../src
meson compile

this might not be needed, build is included in meson config

gcc -fPIC -o test1c test1.c -L../build2 -lrhpix
LD_LIBRARY_PATH=../build2 ./test1c
