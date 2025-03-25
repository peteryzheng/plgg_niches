#!/usr/bin/env bash

pip install anndata scanpy numpy sklearn matplotlib torch datetime
rm -r persist
git clone https://github.com/iancovert/persist.git
pip install -q ./persist
