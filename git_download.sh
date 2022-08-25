#!/bin/bash

from_branch_name="main"
commit_message="fix"

while [ -n "$1" ]; do #setting variables for sub-processes
	case "$1" in
		-from_branch_name) from_branch_name="$2"; shift ;;
		-commit_message) commit_message="$2"; shift ;;
		*) echo "$1 is not an option"; exit 1 ;;
	esac
	shift
done




git reset --hard
git pull origin $from_branch_name