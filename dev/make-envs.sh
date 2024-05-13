#!/bin/bash
echo "Activating base environment"
eval "$(conda shell.bash hook)"
conda activate base

message(){
	let "n = ${#1} + 4"
	echo ""
	for i in $(seq $n); do echo -n "#"; done
	echo ""
	echo -n "# "
	echo -n "$1"
	echo -n " #"
	echo ""
	for i in $(seq $n); do echo -n "#"; done
	echo ""
	echo ""
}

build_myjive(){
	message "CREATING myjive ENVIRONMENT"
	conda env create -f ../ENVIRONMENT.yml
}

build_myjive_dev(){
	message "CREATING myjive-dev ENVIRONMENT"
	conda env create -f ENVIRONMENT-dev.yml -y

	message "ADDING LOCAL PATHS"
	conda activate myjive-dev
	conda develop ~/Storage/git/myjive
	conda deactivate
}

# (re)build myjive environment
if conda env list | grep -q "^myjive "; then
	while true; do
		read -p "myjive environment already exists
Do you want to rebuild it? [Y/n] " yn
		case $yn in
			[Yy]* )
				message "REMOVING myjive ENVIRONMENT"
				conda remove --name myjive --all -y
				build_myjive
				break
				;;
			[Nn]* )
				echo "Skipping myjive build"
				break
				;;
			* )
				echo "Please answer yes or no."
				;;
		esac
	done
else
	build_myjive
fi


# (re)build myjive-dev environment
if conda env list | grep -q "^myjive-dev "; then
	while true; do
		read -p "myjive-dev environment already exists.
Do you want to rebuild it? [Y/n] " yn
		case $yn in
			[Yy]* )
				message "REMOVING myjive-dev ENVIRONMENT"
				conda remove --name myjive-dev --all -y
				build_myjive_dev
				break
				;;
			[Nn]* )
				echo "Skipping myjive-dev build"
				break
				;;
			* )
				echo "Please answer yes or no."
				;;
		esac
	done
else
	build_myjive_dev
fi

