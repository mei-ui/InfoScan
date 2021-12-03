#! /bin/bash
less $1| awk '$3=="transcript"'|grep -Eo 'transcript_id \"\w+\.\w+|transcript_id \"\w+\.\w+\.\w+'|cut -d\" -f 2 |sort|uniq > $2

