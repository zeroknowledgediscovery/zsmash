#!/bin/bash
if [ $# -lt 1 ] ; then
	COMMENT=" updated "
else
	COMMENT=$1
fi

a=$((`awk -F= '{print $2}' version.py  | sed "s/\s*'//g" | awk -F. '{print $NF}'` + 1))
VERSION=`awk -F= '{print $2}' version.py  | sed "s/\s*'//g"`
NEW_VERSION=`awk -F= '{print $2}' version.py  | sed "s/\s*'//g" | awk -F. 'BEGIN{OFS="."}{print $1,$2}'`.$a

echo __version__ = \'$VERSION\'

# git add * -v
# git commit -m "$COMMENT"
# git push
# git tag $VERSION -m "$COMMENT"
# git push --tags

# rm -rf dist
# setup.py needs to be run in manylinux1 docker container
# python3 setup.py sdist
# python3 setup.py bdist_wheel
twine check dist/*
twine upload dist/*

echo __version__ = \'$VERSION\' > version.py
