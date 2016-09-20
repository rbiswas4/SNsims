#!/usr/bin/env bash
pd=$PWD
cd ..
git clone https://github.com/rbiswas4/AnalyzeSN.git
cd AnalyzeSN
python setup.py install --user
cd $pd
