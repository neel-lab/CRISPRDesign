#!/bin/sed
sed -i -e '/^$/d' -e 's/chromosome\_//g' -e 's/\_\(start\|end\)\_/\_/g' $1 
sed -i 's/^>\(.*\)\_\(EN.*$\)/>\2\n\1/g' $1 

