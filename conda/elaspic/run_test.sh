#!/bin/bash

# rsync -av $SRC_DIR/tests ./
cd $SRC_DIR
rm -r tests/__pycache__

py.test -vsx --cache-clear --quick

