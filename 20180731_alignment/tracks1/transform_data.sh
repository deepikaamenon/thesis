#!bin/bash
mkdir original_data
cp -n *W[12]*.txt original_data
matlab < /Users/deepikaa/Desktop/data_desktop/Tracking/DC/wt/Rvs167/2017_04_06_3131/tracks/transform_data.m > /dev/null
