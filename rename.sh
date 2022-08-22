#!/bin/bash

while bjobs | grep -q "$1"; do 
	wait 240
done

echo "done with $1"