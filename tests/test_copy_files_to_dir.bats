#!/usr/bin/env bats

@test "copy files from A to B" {
	#Setup
	ymlPath="$BATS_TEST_DIRNAME/inputs/copy_files.yml"
	touch $ymlPath
	touch ./tests/inputs/small ./tests/inputs/medium ./tests/inputs/large
	echo 'small_word' >> ./tests/inputs/small
	echo 'medium sentence, in one line' >> ./tests/inputs/medium
	echo -e 'large file has multiple\nlines, and can still be copied\nproperly' >> ./tests/inputs/large
	echo 'dirname: "copydir"' >> $ymlPath
	echo 'input_files: [{class: File, path: ./small}, {class: File, path: ./medium}, {class: File, path: ./large}]' >> $ymlPath

	#Testing
	cwltool ./SJARACNe/cwl/copy_files_to_dir.cwl $ymlPath
	run diff ./tests/inputs/small ./copydir/small && diff ./tests/inputs/medium ./copydir/medium && diff ./tests/inputs/large ./copydir/large
	[ "$status" -eq 0 ]
	
	#Cleanup
	rm $ymlPath
	rm ./tests/inputs/small ./tests/inputs/medium ./tests/inputs/large
	rm -r ./copydir
}
