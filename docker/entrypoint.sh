#!/bin/bash
#usermod -u ${LOCAL_USER_ID} imp
#groupmod -g ${LOCAL_GROUP_ID} imp
exec gosu imp "$@"
