#!/bin/sh

homeDIR="$( pwd )"



madgraph="MG5_aMC_v3.3.1.tar.gz"
URL=https://launchpad.net/mg5amcnlo/3.0/3.3.x/+download/$madgraph
echo -n "Install MadGraph (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
	mkdir MG5;
	echo "[installer] getting MadGraph5"; wget $URL 2>/dev/null || curl -O $URL; tar -zxf $madgraph -C MG5 --strip-components 1;
	cd $homeDIR
	rm $madgraph;
fi


