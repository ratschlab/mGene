function [sg_path, envstr] = shogun_settings()

sg_dir=getenv('MGENE_SRC_PATH');
sg_path = deblank(sg_dir);
envstr = sprintf('export LD_LIBRARY_PATH=%s', sg_path);
