#!/bin/bash
USER_ID=${LOCAL_USER_ID:-9001}
usermod -u $USER_ID imp
su imp -c -- $1
