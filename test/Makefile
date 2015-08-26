run: job.sh
	echo "#!/bin/bash" >launch.sh
	echo "# Automatically generated from job.sh - do not edit." >>launch.sh
	echo -n '# @ job_name=' >>launch.sh
	basename ${PWD} >>launch.sh
	cat job.sh >>launch.sh
	llsubmit launch.sh
status:
	llq
#	-u ${USER}
clone: job.sh 
	cp job.sh ../clone
	cp Makefile ../clone
clean:
	rm -f *.[0-9]*.out
	rm -f *.[0-9]*.err
