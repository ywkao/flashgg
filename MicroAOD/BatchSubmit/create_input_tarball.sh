#! /bin/bash

if [ -f package.tar.xz ]; then
    rm ./package.tar.xz
fi

tar -hcJf package.tar.xz --exclude='.git' --exclude='*my*.root' --exclude='test*.root' --exclude='*.tar*' --exclude='BatchSubmit*' --exclude='tasks*' ../../../../../$CMSSW_VERSION