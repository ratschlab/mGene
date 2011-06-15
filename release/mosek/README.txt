To use Mosek with mGene, please follow the these steps:

1. Download the Mosek distribution from http://www.mosek.com ->
Download (tested only with Mosek version 5).  

2. Copy the content of the folders 'mosek/5/tools/platform/*/bin' and
'mosek/5/tools/platform/*/h' in the archive to this directory.

3. Copy the Matlab interface files for Mosek in the archive at
'mosek/5/toolbox/r2007a' into the subdirectory './matlab'.

4. If you want to use Octave with Mosek, go to the subdirectory 'octave'
and run 'make'. It may be necessary to adapt the Makefile (in
particular, depending on the platform, -lmosek may need to be replaced
by -lmosek64).

5. Add the license file (mosek.lic) to this directory. It can be
obtained either by purchasing a mosek license or by using the trial
license part of the download archive in 'mosek/5/licenses/mosek.lic'.

