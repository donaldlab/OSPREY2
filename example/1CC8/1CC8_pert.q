#!/bin/csh -f
# ---------------------------
# our name 
#$ -N test1CC8_pert
#
#
# Operate in current working directory
#$ -cwd
#$ -S /bin/csh
#
#$ -l mem_free=8G
#$ -o 1CC8_pert.out -j y
#
# ---------------------------
set jdep=../../lib
set jars=$jdep/joptimizer.jar:$jdep/architecture-rules-3.0.0-M1.jar:$jdep/colt-1.2.0.jar:$jdep/commons-beanutils-1.6.jar:$jdep/commons-collections-2.1.jar:$jdep/commons-digester-1.6.jar:$jdep/commons-io-1.4.jar:$jdep/commons-lang-2.5.jar:$jdep/commons-logging-1.1.1.jar:$jdep/commons-math3-3.0.jar:$jdep/jdepend-2.9.1.jar:$jdep/junit-3.8.1.jar:$jdep/log4j-1.2.14.jar:$jdep/xml-apis-1.0.b2.jar

java -Xmx4g -cp ../../../build/classes:$jars KStar -c KStar.cfg doDEE System.cfg DEE_pert.cfg
