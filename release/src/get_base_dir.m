function base_dir = get_base_dir()

base_dir = getenv('MGENE_SRC_PATH');
base_dir = deblank(base_dir);
