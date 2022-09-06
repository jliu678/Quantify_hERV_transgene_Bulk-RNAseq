#!/bin/bash

from_branch_name="main"
to_branch_name="main"
commit_message="fix"

while [ -n "$1" ]; do #setting variables for sub-processes
	case "$1" in
		-from_branch_name) from_branch_name="$2"; shift ;;
 		-to_branch_name) to_branch_name="$2"; shift ;;  
		-commit_message) commit_message="$2"; shift ;;
		-f) from_branch_name="$2"; shift ;;
 		-t) to_branch_name="$2"; shift ;;  
		-m) commit_message="$2"; shift ;;
		*) echo "$1 is not an option"; exit 1 ;;
	esac
	shift
done


git add --all
git commit -m "fix"
git push origin $from_branch_name:$to_branch_name