#!/bin/bash
/etc/init.d/postgresql start
while :
do	
	echo "Type search parameters"
	read str
	case "$str" in
		'exit') 
		echo "Exit"
		exit 0
 		;;
	 	*) java -jar vdjdb-1.0-SNAPSHOT.jar -p docker -u docker -d docker $str
		;;
	esac
done
