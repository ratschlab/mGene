Download the most recent version (>0.7.4) of the Shogun toolbox at
http://www.shogun-toolbox.org and compile it for use with Octave or
Matlab (whatever you use). Place the sg.mex* (matlab) or sg.oct (octave)
as well as the libshogun and libshogunui libraries in this directory.

Here is an example of the content of this directory that works with octave:

-rw-r--r--   1 raetsch  fml   239B 28 Mai 14:36 README.txt
-rwxr-xr-x   1 raetsch  fml   2,3M 30 Mai 10:24 libshogun.4.0.dylib*
lrwxr-xr-x   1 raetsch  fml    19B 30 Mai 10:24 libshogun.4.dylib@ -> libshogun.4.0.dylib
-rw-r--r--   1 raetsch  fml   9,8M 30 Mai 10:24 libshogun.a
lrwxr-xr-x   1 raetsch  fml    19B 30 Mai 10:24 libshogun.dylib@ -> libshogun.4.0.dylib
-rwxr-xr-x   1 raetsch  fml   739K 17 Mai 17:18 libshogunui.1.2.dylib*
lrwxr-xr-x   1 raetsch  fml    21B 17 Mai 17:18 libshogunui.1.dylib@ -> libshogunui.1.2.dylib
-rwxr-xr-x   1 raetsch  fml   747K 30 Mai 10:24 libshogunui.2.0.dylib*
lrwxr-xr-x   1 raetsch  fml    21B 30 Mai 10:24 libshogunui.2.dylib@ -> libshogunui.2.0.dylib
-rw-r--r--   1 raetsch  fml   2,8M 30 Mai 08:58 libshogunui.a
lrwxr-xr-x   1 raetsch  fml    21B 30 Mai 10:24 libshogunui.dylib@ -> libshogunui.2.0.dylib
-rwxr-xr-x   1 raetsch  fml   183K 30 Mai 10:30 sg.oct*
