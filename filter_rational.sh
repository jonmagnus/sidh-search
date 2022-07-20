#!/bin/bash
grep -E -e "(i\*0.*){2}" -e "[{}]" $1 | neato -Tpdf -Goverlap=scale -o $2
