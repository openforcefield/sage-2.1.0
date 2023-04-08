#!/bin/bash -x

mkdir -p /tmp/dir/fb_193
cd /tmp/dir
while [ -f state.COPYING ] ; do sleep 100 ; done
if [ ! -f fb_193.tar.gz ] ; then 
	touch state.COPYING
	cp /dir/dir/dir/conda-env/fb_193.tar.gz .
    tar xzf fb_193.tar.gz -C fb_193
    rm state.COPYING
fi
source fb_193/bin/activate
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

work_queue_worker $@
