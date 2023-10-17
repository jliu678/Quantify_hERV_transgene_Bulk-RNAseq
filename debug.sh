. ~/scratch/erv/herv_project/timed.sh

get_r2name(){
	local r2id="r2"
	if echo $1 | grep -q "R1"; then r2id="R2"; fi
	echo $1 | sed "s/r1/$r2id/i"
}

check_name(){
	#echo $1 | grep -iq "r1"
	#local is_r1=$?
    echo $1
  if echo $1 | grep -iq "r1"
  then
    local is_r1=true
  else
    local is_r1=false
  fi
  timed_print "===!!! is_r1 is $is_r1==="
  
  if [ -e "$(get_r2name $1)" ]
  then
    local r2_exist=true
  else
    local r2_exist=false
  fi
  timed_print "===!!! $(get_r2name $1) is $r2_exist==="

  if $is_r1 && $r2_exist
  then
    timed_print "===!!! return is 0===" 
    return 0
  else
    timed_print "===!!! return is 1===" 
    return 1
  fi
}

check_name ../raw_data/caki127_3-k25-both.R1.fq.gz