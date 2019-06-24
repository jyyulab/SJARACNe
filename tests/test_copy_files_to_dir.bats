#!/usr/bin/env bats

@test "copy files from A to B" {
	#Setup
	ymlPath="$BATS_TEST_DIRNAME/inputs/copy_files.yml"
	touch $ymlPath
	touch ./tests/inputs/small ./tests/inputs/medium ./tests/inputs/large
	echo 'dirname: "test"' >> $ymlPath
	echo 'input_files: [{class: File, path: ./small}, {class: File, path: ./medium}, {class: File, path: ./large}]' >> $ymlPath

	#Testing
	cwltool ./SJARACNe/cwl/copy_files_to_dir.cwl $ymlPath | tee results
	
	#Cleanup
	rm $ymlPath
	rm ./tests/inputs/small ./tests/inputs/medium ./tests/inputs/large
	rm -r test
}
