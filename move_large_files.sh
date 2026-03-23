SOURCE_DIR="./"
DEST_DIR="../blog_03_12_2025_largefiles_huggingface"
MIN_SIZE="+5M"
# Ensure the destination directory exists
#mkdir -p "$DEST_DIR"
#
#find "$SOURCE_DIR" -type f -size "$MIN_SIZE" \
    #\( -name '*.xtc' -o \
       #-name '*.dcd' -o \
       #-name '*.nc' -o \
       #-name '*.pdb' -o \
       #-name '*.densities' -o \
       #-name '*.gbw' -o \
       #-name '*.molden' -o \
       #-name '*.model' -o \
       #-name '*.pt' \
    #\) \
    #-exec cp --parents {} "$DEST_DIR" \;
#
find "$SOURCE_DIR" -type f -name '*.pdb' -exec cp --parents {} "$DEST_DIR" \;
