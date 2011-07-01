
basedir="/home/galaxy/svn/projects/mGene_core/"

envstr=`octave --eval "cd $basedir; [extra_path shogun_envstr] = shogun_settings(); fprintf(1, '%s', shogun_envstr);" 2>/dev/null | grep LD_LIBRARY_PATH`

echo $envstr
