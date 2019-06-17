#!/usr/bin/env bash

./SJARACNe/sjaracne.py lsf -j ./SJARACNe/config/config_cwlexec.json -e /research/rgs01/project_space/yu3grp/Network_JY/yu3grp/LiangDing/chemoresistance-brca/TCGA/TCGA_V201503_TN_139/project_2019-06-10/SJAR/project_2019-06-10/input.exp -g /research/rgs01/project_space/yu3grp/Network_JY/yu3grp/LiangDing/chemoresistance-brca/Metabric/TN_271/project_2019-04-19/SJAR/project_2019-04-19/qingfei_lists/tf.txt -n 2 -o ./test_data/outputs/cwl/cwltool/SJARACNE_out.final

