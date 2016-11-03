#!/usr/bin/env bash
pd=$PWD
cd ..
git clone https://github.com/rbiswas4/AnalyzeSN.git
cd AnalyzeSN
python setup.py install --user
# source ./continuous_integration/setup_env.sh
# pip install -r ./continuous_integration/pip_requirements.txt
#python setup.py install --user
cd $pd
