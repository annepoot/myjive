# Check existence of all files
for number in {1..9..2}; do
	FILE=hole-$number.msh
	if test -f "$FILE"; then
		echo "$FILE exists."
	else
		echo -e "$FILE is missing.\nScript execution is terminated."
		exit 1
	fi
done

echo -e ""

# Check existence of all folders, or create them
for number in 4 9; do
	FOLDER="q$number"
	if test -d "$FOLDER"; then
		echo "$FOLDER already exists."
	else	
		echo -e "$FOLDER does not yet exist\nCreating $FOLDER folder."
		mkdir "$FOLDER"
	fi
done

echo ""

# Create the q4 meshes
for number in {1..9..2}; do
	FILE=hole-$number.msh
	FOLDER=q4/
	cp "$FILE" "$FOLDER"

	OFILE=hole-$number-r0.msh
	cp "$FOLDER$FILE" "$FOLDER$OFILE"
	
	echo ""
	for rlevel in {1..4}; do
		OFILE=hole-$number-r$((rlevel-1)).msh
		RFILE=hole-$number-r$rlevel.msh

		echo -e "Refining the mesh from $OFILE to $RFILE"
		gmsh -refine "$FOLDER$OFILE" -format msh22 -order 1 -v 3 -o "$FOLDER$RFILE"
	done
done

echo ""

# Create the q9 meshes
for number in {1..9..2}; do
	OFOLDER=q4/
	RFOLDER=q9/

	echo ""

	FILE=hole-$number.msh
	echo -e "Changing the element order of $FILE from q4 to q9"
	gmsh -2 "$OFOLDER$FILE" -format msh22 -order 2 -v 3 -o "$RFOLDER$FILE"

	for rlevel in {0..4}; do
		FILE=hole-$number-r$rlevel.msh
		echo -e "Changing the element order of $FILE from q4 to q9"
		gmsh -2 "$OFOLDER$FILE" -format msh22 -order 2 -v 3 -o "$RFOLDER$FILE"
	done
done

echo -e "\n\nMesh generation procedure finished."
