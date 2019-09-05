#!/bin/bash
#change the number after  -j change the number of files to be processed.
parallel --verbose -j 2 prefetch {} ::: $(cut -f1 remain_list ) >>sra_download_remain.log
wait
exit
