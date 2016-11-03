#!/usr/bin/env bash
pd=$PWD
cd ..
git clone https://github.com/esheldon/fitsio.git
cd fitsio
pyton setup.py install --user
cd ..
git clone https://github.com/rbiswas4/OpSimSummary.git
cd OpSimSummary
python setup.py install --user
cd $pd
