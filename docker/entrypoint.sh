#!/bin/bash

USER_ID=${LOCAL_USER_ID:1000}
GROUP_ID=${LOCAL_GROUP_ID:1000}
usermod -u $USER_ID imp
groupmod -g $GROUP_ID imp
exec gosu imp "$@"
